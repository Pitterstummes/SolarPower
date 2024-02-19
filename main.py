import numpy as np
import matplotlib.pyplot as plt
import datetime
import math
import ephem


def calculate_sun_position(latitude, longitude, date):
    # Calculate the sun's position using the given coordinates and date
    observer = ephem.Observer()
    observer.lat = str(latitude)
    observer.lon = str(longitude)
    observer.date = date

    sun = ephem.Sun()
    sun.compute(observer)

    azimuth = np.rad2deg(sun.az)
    elevation = np.rad2deg(sun.alt)
    return azimuth, elevation


def calculate_delta_angle(sun_azimut, panel_orientation, sun_elevation, panel_tilt):
    # Calculate the angle between the sun and the solar panel
    # Convert to radians
    sun_azimut = np.deg2rad(sun_azimut)
    panel_orientation = np.deg2rad(panel_orientation)
    sun_elevation = np.deg2rad(sun_elevation)
    panel_tilt = np.deg2rad(panel_tilt)
    # Convert to spherical coordinates, r = 1
    theta_sun = np.pi / 2 - sun_elevation
    theta_panel = np.pi / 2 - panel_tilt
    phi_sun = sun_azimut
    phi_panel = panel_orientation
    # Calculate dot product of the two vectors
    dot_product = np.sin(theta_sun) * np.sin(theta_panel) * np.cos(
        phi_sun - phi_panel
    ) + np.cos(theta_sun) * np.cos(theta_panel)
    # Calculate angle in degrees
    angle = np.rad2deg(np.arccos(dot_product))
    return angle


def create_timestamps(date, time_step):
    # Create an array of timestamps for one day
    start_time = datetime.datetime(date.year, date.month, date.day)
    end_time = start_time + datetime.timedelta(days=1)
    num_steps = int((end_time - start_time) / time_step)
    timestamps = np.array([start_time + i * time_step for i in range(num_steps)])
    return timestamps


def get_sun_position(
    latitude, longitude, date, time_step, positive_elevation_mask=True
):
    # Calculate the sun's position for each timestamp
    timestamps = create_timestamps(date, time_step)
    azimuths = np.empty(len(timestamps))
    elevations = np.empty(len(timestamps))
    for i, timestamp in enumerate(timestamps):
        azimuth, elevation = calculate_sun_position(latitude, longitude, timestamp)
        azimuths[i] = azimuth
        elevations[i] = elevation
    # Remove negative elevation values
    if positive_elevation_mask:
        positive_elevation_mask = np.array(elevations) >= 0
        azimuths = np.array(azimuths)[positive_elevation_mask]
        elevations = np.array(elevations)[positive_elevation_mask]
        timestamps = np.array(timestamps)[positive_elevation_mask]
    return azimuths, elevations, timestamps


def plot_sun_position(azimuths, elevations, timestamps):
    # Plot the sun's position
    fig, ax = plt.subplots()

    def time_to_decimal(time):
        return time.hour + time.minute / 60

    ax.scatter(
        azimuths,
        elevations,
        c=np.linspace(
            time_to_decimal(timestamps[0]),
            time_to_decimal(timestamps[-1]),
            len(timestamps),
        ),
        cmap="hsv",
        label="Azimuth vs Elevation",
    )
    ax.set_xlabel("Azimuth in degrees")
    ax.set_ylabel("Elevation in degrees")
    titletext = (
        "Sun Position for Lat: "
        + str(latitude)
        + " and Lon: "
        + str(longitude)
        + " on "
        + str(date.date())
    )
    ax.set_title(titletext)

    def azimuth_to_cardinal(azimuth):
        return azimuth / 45

    def cardinal_to_azimuth(cardinal):
        return cardinal * 45

    secax_x = ax.secondary_xaxis(
        "top", functions=(azimuth_to_cardinal, cardinal_to_azimuth)
    )

    fig.colorbar(
        mappable=ax.collections[0],
        format=lambda t, pos: f"{math.floor(t):02.0f}:{(t*60)%60:02.0f}",
    )

    plt.grid(True)
    plt.show()


def get_pv_angle(azimuths, elevations, panel_azimuth, panel_elevation):
    # Calculate the angle between the sun and the solar panel for each timestamp
    angles = np.empty((len(panel_azimuth), len(azimuths)))
    for num_panel in range(len(panel_azimuth)):
        for i in range(len(azimuths)):
            angles[num_panel, i] = calculate_delta_angle(
                azimuths[i],
                panel_azimuth[num_panel],
                elevations[i],
                panel_elevation[num_panel],
            )
    return angles, azimuths


def calc_power_ratio(angles):
    # Calculate the power ratio for each angle
    power_ratio = np.zeros(len(angles[0][0]))
    for i in range(len(power_ratio)):
        for j in range(len(angles)):
            if angles[0][j][i] < 90:
                power_ratio[i] += np.cos(np.deg2rad(angles[0][j][i]))
    return power_ratio


def calc_irradiance(elevations):
    # Calculate the irradiance for each elevation
    irradiance = np.zeros(len(elevations))
    for i in range(len(elevations)):
        irradiance[i] = np.sin(np.deg2rad(elevations[i]))
    return irradiance


def plot_panel_angle(
    azimuths,
    angles,
    elevations,
    timestamps,
    panel_azimuth,
    panel_elevation,
):
    # Plot the angle between the sun and the solar panel
    fig, ax = plt.subplots()
    for i in range(len(panel_azimuth)):
        x = []
        y = []
        for j in range(len(azimuths)):
            if angles[0][i, j] < 90:
                x.append(azimuths[j])
                y.append(angles[0][i][j])
        ax.scatter(
            x,
            y,
            label="Panel " + str(i + 1),
        )
    power_ratio = calc_power_ratio(angles)
    irradiance = calc_irradiance(elevations)
    irrad_power_ratio = power_ratio * irradiance
    irrad_power_ratio_norm = irrad_power_ratio / np.max(irrad_power_ratio) * 90
    ax.scatter(azimuths, irrad_power_ratio_norm, label="Power Ratio")
    ax.set_xlabel("Azimuth in degrees")
    ax.set_ylabel("Angle in degrees")
    titletext = "Angle between Sun and Solar Panel"
    ax.set_title(titletext)
    ax.legend()
    plt.grid(True)
    plt.show()


# Coordinates & date
latitude = 48.8
longitude = 9
date = datetime.datetime.now() + datetime.timedelta(days=0)  # Add days
time_step = datetime.timedelta(minutes=1)
azimuths_sun, elevations_sun, timestamps = get_sun_position(
    latitude, longitude, date, time_step, True
)
plot_sun_position(azimuths_sun, elevations_sun, timestamps)
panel_azimuths = np.array([96, 96 + 180])
panel_elevations = np.array([30, 45])
angles = get_pv_angle(azimuths_sun, elevations_sun, panel_azimuths, panel_elevations)
plot_panel_angle(
    azimuths_sun, angles, elevations_sun, timestamps, panel_azimuths, panel_elevations
)
