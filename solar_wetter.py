"""PV yield with real weather: SMA data, Open-Meteo weather, pvlib model, calibration, forecast.

    from solar_wetter import *
    tage = sma_tage()                       # cleaned daily yields (with flags)
    wetter = wetter_archiv("2023-07-20", "2025-12-31")
    stunden = modell(wetter)                # hourly kWh per array (before calibration)
    k, vergleich = kalibrieren(tage, stunden)
    prognose(k)                             # next 7 days, hourly kWh

Model: sun position and irradiance on the two tilted arrays (pvlib, Hay-Davies),
direct sunlight blocked by the measured horizon (terrain_example_1.csv), cell
temperature (Faiman) with -0.4 %/K. All remaining losses (inverter, cables, dirt,
shading by snow, ...) are one calibration factor k, fitted to the measured days.

Weather: Open-Meteo (free, no key). Archive = reanalysis (ERA5 + others) for the past,
forecast API for the next days. Only the rounded coordinates are sent, no address.
"""

import os

import numpy as np
import pandas as pd
import pvlib
import requests

HERE = os.path.dirname(os.path.abspath(__file__))
DATEN = os.path.join(HERE, "Daten")
SMA = os.path.join(DATEN, "SMA Export")

# Location and arrays. Azimuths fitted to the 5-minute data Jul-Dec 2025 give 90-100 deg (east)
# and 275-285 deg (west), so the values from main.ipynb are kept.
LAT, LON = 48.80, 9.01          # Leonberg, rounded
ARRAYS = {                       # name: (azimuth deg, 0=N 90=E; tilt deg; kWp)
    "Ost": (96, 30, 3.32),       # 8 modules x 415 Wp; tilt about 30 deg (Paul, 09.10.2026)
    "West": (276, 45, 4.98),     # 12 modules x 415 Wp; tilt 45 deg (Paul, 09.10.2026)
}
GAMMA = -0.004                   # power temperature coefficient per K
HORIZON_FILE = os.path.join(HERE, "terrain_example_1.csv")


# ------------------------------------------------------------------ measured data
def sma_tage():
    """Daily SMA values 2023-.. with flags. Corrections:
    - days with generation 0 AND feed-in 0 in 2023: system switched off (not registered
      yet) -> flag 'aus', excluded from calibration
    - logger outage: generation 0 but feed-in > 0 -> the day was added to the next day.
      The next day is split: its own generation = own feed-in + typical direct use,
      the rest belongs to the outage day. Both flagged 'korrigiert'."""
    files = sorted(f for f in os.listdir(SMA) if f.endswith("_Energy.csv"))
    d = pd.concat(pd.read_csv(os.path.join(SMA, f), sep=";", skiprows=1, decimal=",", na_values="-") for f in files)
    d["Zeitraum"] = pd.to_datetime(d["Zeitraum"], dayfirst=True)
    d = d.set_index("Zeitraum").sort_index()
    d.columns = ["erzeugung", "verbrauch", "direkt", "einspeisung", "bezug", "bat_laden", "bat_entladen"]
    d["flag"] = "ok"
    first = d.index[0]
    d.loc[first, "flag"] = "aus"                                   # commissioning day
    off = (d.erzeugung.fillna(0) == 0) & (d.einspeisung.fillna(0) == 0)
    d.loc[off, "flag"] = "aus"
    outage = (d.erzeugung.fillna(0) == 0) & (d.einspeisung.fillna(0) > 0)
    for day in d.index[outage]:
        nxt = day + pd.Timedelta(days=1)
        if nxt not in d.index:
            continue
        neighbours = d.loc[day - pd.Timedelta(days=5):nxt + pd.Timedelta(days=5)]
        typical_direct = neighbours.loc[neighbours.flag == "ok", "direkt"].drop([nxt], errors="ignore").median()
        total = d.loc[nxt, "erzeugung"]
        own = d.loc[nxt, "einspeisung"] + typical_direct
        d.loc[nxt, "erzeugung"] = own
        d.loc[day, "erzeugung"] = total - own
        d.loc[[day, nxt], "flag"] = "korrigiert"
    return d


# ------------------------------------------------------------------ weather
def _openmeteo(url, params):
    p = dict(latitude=LAT, longitude=LON, timezone="UTC",
             hourly="shortwave_radiation,direct_normal_irradiance,diffuse_radiation,temperature_2m,wind_speed_10m")
    p.update(params)
    r = requests.get(url, params=p, timeout=60)
    r.raise_for_status()
    h = r.json()["hourly"]
    w = pd.DataFrame(h).rename(columns=dict(shortwave_radiation="ghi", direct_normal_irradiance="dni",
                                            diffuse_radiation="dhi", temperature_2m="temp", wind_speed_10m="wind"))
    w["time"] = pd.to_datetime(w["time"], utc=True)
    w = w.set_index("time")
    w["wind"] = w["wind"] / 3.6          # km/h -> m/s
    return w.astype(float)


def wetter_archiv(start, end, cache=True):
    """Hourly weather of the past (values are means of the PRECEDING hour)."""
    path = os.path.join(DATEN, f"wetter_{start}_{end}.csv")
    if cache and os.path.exists(path):
        return pd.read_csv(path, index_col=0, parse_dates=True)
    w = _openmeteo("https://archive-api.open-meteo.com/v1/archive", dict(start_date=start, end_date=end))
    os.makedirs(DATEN, exist_ok=True)
    w.to_csv(path)
    return w


def wetter_prognose(days=7):
    return _openmeteo("https://api.open-meteo.com/v1/forecast", dict(forecast_days=days))


# ------------------------------------------------------------------ model
def _horizon():
    t = np.genfromtxt(HORIZON_FILE, delimiter=",")
    az = np.concatenate([[0], t[:, 0], [360]])
    el = np.concatenate([[0], t[:, 1], [0]])
    order = np.argsort(az)
    return az[order], el[order]


def modell(w):
    """Hourly DC energy (kWh) per array before calibration, index = hour (UTC, end of interval)."""
    mid = w.index - pd.Timedelta(minutes=30)                     # sun position in the middle of the hour
    sp = pvlib.solarposition.get_solarposition(mid, LAT, LON)
    dni_extra = pvlib.irradiance.get_extra_radiation(mid)
    haz, hel = _horizon()
    blocked = sp["apparent_elevation"].values < np.interp(sp["azimuth"].values % 360, haz, hel)
    dni = np.where(blocked, 0.0, w["dni"].values)
    out = pd.DataFrame(index=w.index)
    for name, (az, tilt, kwp) in ARRAYS.items():
        poa = pvlib.irradiance.get_total_irradiance(tilt, az, sp["apparent_zenith"].values, sp["azimuth"].values,
                                                    dni, w["ghi"].values, w["dhi"].values,
                                                    dni_extra=dni_extra.values, model="haydavies")
        g = np.nan_to_num(np.asarray(poa["poa_global"], dtype=float)).clip(min=0)
        tcell = pvlib.temperature.faiman(g, w["temp"].values, w["wind"].values)
        out[name] = g / 1000 * kwp * (1 + GAMMA * (tcell - 25))
    out["gesamt"] = out[list(ARRAYS)].sum(axis=1)
    return out


def tageswerte(stunden):
    local = stunden.tz_convert("Europe/Berlin")
    return local.groupby(local.index.date).sum().rename_axis("tag").set_index(pd.DatetimeIndex(
        pd.to_datetime(sorted(set(local.index.date)))))


# ------------------------------------------------------------------ calibration
def kalibrieren(tage, stunden, monatlich=True):
    """Calibration factor(s) by least squares over valid days. Winter has extra losses
    (low sun, snow, unmeasured horizon), so by default one factor per month.
    Returns k (Series month -> factor) and a comparison table."""
    mod = tageswerte(stunden)["gesamt"]
    df = pd.DataFrame({"gemessen": tage["erzeugung"], "flag": tage["flag"], "modell_roh": mod}).dropna(subset=["modell_roh"])
    ok = df[(df.flag == "ok") & df.gemessen.notna()]
    k_all = float((ok.gemessen * ok.modell_roh).sum() / (ok.modell_roh ** 2).sum())
    if monatlich:
        k = ok.groupby(ok.index.month).apply(lambda s: (s.gemessen * s.modell_roh).sum() / (s.modell_roh ** 2).sum())
        k = k.reindex(range(1, 13)).fillna(k_all)
    else:
        k = pd.Series(k_all, index=range(1, 13))
    df["modell"] = df["modell_roh"] * df.index.month.map(k).values
    return k, df

def guete(df):
    ok = df[(df.flag == "ok") & df.gemessen.notna()]
    err = ok.modell - ok.gemessen
    r2 = 1 - (err ** 2).sum() / ((ok.gemessen - ok.gemessen.mean()) ** 2).sum()
    return dict(tage=len(ok), r2=r2, mae_kwh=err.abs().mean(), mae_prozent=err.abs().sum() / ok.gemessen.sum() * 100,
                bias_prozent=err.sum() / ok.gemessen.sum() * 100)


def prognose(k, days=7):
    """Hourly forecast (kWh, local time) for the next days; k = factor or Series month -> factor."""
    st = modell(wetter_prognose(days)).tz_convert("Europe/Berlin")
    factor = st.index.month.map(k).values if isinstance(k, pd.Series) else k
    return st.mul(factor, axis=0)
