import numpy as np
import matplotlib.pyplot as plt
import datetime
import math
import ephem


def calculate_sun_position(latitude, longitude, date):
    # Calculate the sun's position using the given coordinates and date
    # You can use any algorithm or library to calculate the sun's position
    # Here's a simple example using the ephem library:

    observer = ephem.Observer()
    observer.lat = str(latitude)
    observer.lon = str(longitude)
    observer.date = date

    sun = ephem.Sun()
    sun.compute(observer)

    azimuth = math.degrees(sun.az)
    elevation = math.degrees(sun.alt)

    return azimuth, elevation


# Example usage
latitude = 48.8  # Latitude of the coordinates
longitude = 9  # Longitude of the coordinates
date = datetime.datetime.now() + datetime.timedelta(days=0)  # Current date and time


# Set the time step to 15 minutes
time_step = datetime.timedelta(minutes=15)

# Create a list of timestamps for each 15 minute step of the day
timestamps = []
start_time = datetime.datetime(date.year, date.month, date.day)
end_time = start_time + datetime.timedelta(days=1)
current_time = start_time
while current_time < end_time:
    timestamps.append(current_time)
    current_time += time_step

# Calculate the sun's position for each timestamp
azimuths = []
elevations = []
for timestamp in timestamps:
    azimuth, elevation = calculate_sun_position(latitude, longitude, timestamp)
    azimuths.append(azimuth)
    elevations.append(elevation)

x = azimuths
y = elevations
if True:
    positive_elevation_mask = np.array(elevations) > -5
    x = np.array(x)[positive_elevation_mask]
    y = np.array(y)[positive_elevation_mask]
    timestamps = np.array(timestamps)[positive_elevation_mask]

# Plot
plot = True
if plot:
    fig, ax = plt.subplots(layout="constrained")

    def time_to_decimal(time):
        return time.hour + time.minute / 60

    ax.scatter(
        x,
        y,
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


def calculate_angle(azimuth1, elevation1, azimuth2, elevation2):
    # Calculate the angle between two points in azimuth-elevation coordinates
    # Calculate spherical coordinates
    phi1 = math.radians(azimuth1)
    theta1 = math.radians(90 - elevation1)
    phi2 = math.radians(azimuth2)
    theta2 = math.radians(90 - elevation2)

    angle_rad = math.acos(
        math.sin(phi1) * math.sin(phi2) * math.cos(theta1 - theta2)
        + math.cos(phi1) * math.cos(phi2)
    )

    # Convert the angle to degrees
    angle_deg = math.degrees(angle_rad)

    return angle_deg


pv_azimuth_1_deg = 96.0
pv_elevation_1_deg = 30.0
pv_azimuth_2_deg = 96.0 + 180
pv_elevation_2_deg = 45.0

anglediff1 = []
anglediff2 = []
total_power_ratio = []
for i in range(len(x)):
    angle1 = calculate_angle(x[i], y[i], pv_azimuth_1_deg, pv_elevation_1_deg)
    angle2 = calculate_angle(x[i], y[i], pv_azimuth_2_deg, pv_elevation_2_deg)
    anglediff1.append(angle1)
    anglediff2.append(angle2)
    if angle1 < 90 and angle2 < 90:
        total_power_ratio.append(
            (math.cos(math.radians(angle1)) + math.cos(math.radians(angle2))) * 100
        )
    elif angle1 < 90:
        total_power_ratio.append(math.cos(math.radians(angle1)) * 100)
    elif angle2 < 90:
        total_power_ratio.append(math.cos(math.radians(angle2)) * 100)
    else:
        total_power_ratio.append(0)
plt.scatter(x, anglediff1, label="PV 1: East")
plt.scatter(x, anglediff2, label="PV 2: West")
plt.scatter(x, total_power_ratio, label="Total Power Ratio (a.u.)")
plt.xlabel("Azimuth in degrees")
plt.ylabel("Angle in degrees")
plt.title("Angle between Sun and PV")
plt.grid(True)
plt.legend()
plt.show()
