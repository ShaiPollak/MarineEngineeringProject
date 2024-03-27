import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import pandas as pd
import folium

def plot_ssh(H, start_pixel_index, finish_line_index, pix_dist):
    x_min, x_max = 0, pix_dist*start_pixel_index
    y_min, y_max = 0, pix_dist * (finish_line_index-start_line_index)

    plt.figure(figsize=(10, 6))
    plt.imshow(H, cmap='jet', aspect='auto', extent=[x_min, x_max, y_max, y_min])
    plt.colorbar(label='SSH (m)')
    plt.xlabel('(m) - East')
    plt.ylabel('(m) - North')
    plt.title('Sea Surface Height')
    plt.grid(True)
    plt.show()

# Path to the NetCDF file
nc_file_path = 'Data/SWOT_L2_LR_SSH_Unsmoothed_555_015_20230617T211922_20230617T220935_PIB0_01.nc'

start_line_index = 21930
finish_line_index = 22000
start_pixel_index = 1
finish_pixel_index = 240  # Assuming this is the correct variable name instead of start_pixel_index twice


# Open the NetCDF file
with Dataset(nc_file_path, "r") as netcdf_file:
    # Extracting data
    ssh_left = netcdf_file.groups['left'].variables['ssh_karin_2'][start_line_index:finish_line_index, :]
    ssh_right = netcdf_file.groups['right'].variables['ssh_karin_2'][start_line_index:finish_line_index, :]

    latitude_value_left = netcdf_file.groups['left'].variables['latitude'][start_line_index:finish_line_index]
    latitude_value_right = netcdf_file.groups['right'].variables['latitude'][start_line_index:finish_line_index]

    longitude_value_left = netcdf_file.groups['left'].variables['longitude'][start_line_index:finish_line_index]
    longitude_value_right = netcdf_file.groups['right'].variables['longitude'][start_line_index:finish_line_index]

        # Flattening SSH data
    ssh_left_flat = ssh_left.flatten()
    ssh_right_flat = ssh_right.flatten()

        # Repeat or tile latitude and longitude to match the flattened SSH data
    latitude_expanded_left = np.repeat(latitude_value_left, ssh_left.shape[1])
    longitude_expanded_left = np.tile(longitude_value_left, len(latitude_value_left))

    latitude_expanded_right = np.repeat(latitude_value_right, ssh_right.shape[1])
    longitude_expanded_right = np.tile(longitude_value_right, len(latitude_value_right))

    print("Shape of ssh_left_flat:", ssh_left_flat.shape)
    print("Shape of longitude_expanded_left:", longitude_expanded_left.shape)
    print("Shape of latitude_expanded_left:", latitude_expanded_left.shape)

    print("Shape of ssh_right_flat:", ssh_right_flat.shape)
    print("Shape of longitude_expanded_right:", longitude_expanded_right.shape)
    print("Shape of latitude_expanded_right:", latitude_expanded_right.shape)

    # Create DataFrames
    data_left = pd.DataFrame({
        'Latitude': latitude_expanded_left,
        'Longitude': longitude_expanded_left,
        'SSH': ssh_left_flat
    })

    data_right = pd.DataFrame({
        'Latitude': latitude_expanded_right,
        'Longitude': longitude_expanded_right,
        'SSH': ssh_right_flat
    })

    # Combine left and right data
    data = pd.concat([data_left, data_right])

    # Normalize SSH values for coloring
    data['SSH_normalized'] = (data['SSH'] - data['SSH'].min()) / (data['SSH'].max() - data['SSH'].min())

    # Create a folium map centered around the mean latitude and longitude
    map = folium.Map(location=[data['Latitude'].mean(), data['Longitude'].mean()], zoom_start=5)

    # Add points to the map
    for _, row in data.iterrows():
        folium.CircleMarker(
            location=[row['Latitude'], row['Longitude']],
            radius=5,
            fill=True,
            color=plt.cm.viridis(row['SSH_normalized']),
            fill_color=plt.cm.viridis(row['SSH_normalized']),
            fill_opacity=0.7
        ).add_to(map)

    # Save the map to an HTML file
    map.save('ssh_map.html')

