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


def calculate_delta_angle(azimut_sun, azimut_panel, elevation_sun, elevation_panel):
    # Calculate the angle between the sun and the solar panel
    # Handle degree or radian input
    if azimut_sun > 2 * np.pi:
        azimut_sun = np.deg2rad(azimut_sun)
    if azimut_panel > 2 * np.pi:
        azimut_panel = np.deg2rad(azimut_panel)
    if elevation_sun > 2 * np.pi:
        elevation_sun = np.deg2rad(elevation_sun)
    if elevation_panel > 2 * np.pi:
        elevation_panel = np.deg2rad(elevation_panel)
    # Convert to spherical coordinates, r = 1
    theta_sun = np.pi / 2 - elevation_sun
    print("theta_sun", theta_sun, "elevation_sun", elevation_sun)
    theta_panel = np.pi / 2 - elevation_panel
    phi_sun = azimut_sun
    phi_panel = azimut_panel
    # Calculate dot product of the two vectors
    dot_product = np.sin(phi_sun) * np.sin(phi_panel) * np.cos(
        theta_sun - theta_panel
    ) + np.cos(phi_sun) * np.cos(phi_panel)
    # Calculate angle
    angle = np.rad2deg(np.arccos(dot_product))

    return angle


# Coordinates & date
latitude = 48.8
longitude = 9
date = datetime.datetime.now() + datetime.timedelta(days=0)  # Add days

print(calculate_delta_angle(30, 30, 0, -20))
