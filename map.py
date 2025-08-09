import folium
from folium import features
import matplotlib.pyplot as plt
import numpy as np
import io
import base64

import folium
import pandas as pd

import plot_ssh

# Assuming you have the latitude, longitude, and SSH values in the variables
# latitude_value_left, longitude_value_left, and ssh_left (similarly for _right),
# let's prepare the data:
data_left = pd.DataFrame({
    'Latitude': plot_ssh.latitude_value_left,
    'Longitude': plot_ssh.longitude_value_left,
    'SSH': plot_ssh.latitude_value_left.ssh_left,
})

data_right = pd.DataFrame({
    'Latitude': plot_ssh.latitude_value_right,
    'Longitude': plot_ssh.longitude_value_right,
    'SSH': plot_ssh.latitude_value_left.ssh_right,
})

# Combine the left and right data
data = pd.concat([data_left, data_right])

# Normalize SSH values for coloring (optional, depending on your SSH data range)
data['SSH_normalized'] = (data['SSH'] - data['SSH'].min()) / (data['SSH'].max() - data['SSH'].min())

# Create a folium map centered around the mean latitude and longitude
map = folium.Map(location=[data['Latitude'].mean(), data['Longitude'].mean()], zoom_start=5)

# Add points to the map
for _, row in data.iterrows():
    folium.CircleMarker(
        location=[row['Latitude'], row['Longitude']],
        radius=5,  # Adjust size of the circle markers
        fill=True,
        color=plt.cm.viridis(row['SSH_normalized']),  # Use matplotlib's colormap
        fill_color=plt.cm.viridis(row['SSH_normalized']),
        fill_opacity=0.7
    ).add_to(map)

# Save the map to an HTML file
map.save('ssh_map.html')