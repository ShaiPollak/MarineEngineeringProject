import numpy as np
import folium
from netCDF4 import Dataset
import branca.colormap as cm
import matplotlib.pyplot as plt


def process_ssh_chunk(ssh_data, lat_values, lon_values, map_obj):

    """
    Process a chunk of SSH data and add it to the folium map, handling NaN values explicitly
    and ensuring not to perform operations on empty arrays.
    """
    # Filter out NaN values from latitude or longitude
    valid_indices = (~np.isnan(lat_values)) & (~np.isnan(lon_values))

    # Use valid_indices to filter out NaN values
    valid_lat = lat_values[valid_indices]
    valid_lon = lon_values[valid_indices]
    valid_ssh = ssh_data[valid_indices]

    # Check if the resulting valid_ssh array is empty
    if valid_ssh.size == 0:
        print("Filtered SSH data is empty for this chunk. Skipping...")
        return

    # Normalize the valid SSH values for coloring
    ssh_normalized = (valid_ssh - np.nanmin(valid_ssh)) / (np.nanmax(valid_ssh) - np.nanmin(valid_ssh))

    #Define a Con. Color Map
    ssh_range_min = np.nanmin(valid_ssh)
    ssh_range_max = np.nanmax(valid_ssh)
    color_map = cm.linear.YlGnBu_09.scale(ssh_range_min, ssh_range_max)
    color_map.caption = 'Sea Surface Height (m)'  # Add a caption for the legend
    map_obj.add_child(color_map) #Add the Colormap to the Map as a Legend

    #Color and add each dot to the map
    for lat, lon, ssh_norm in zip(valid_lat, valid_lon, ssh_normalized):
        # Skip plotting if any value is NaN - this check is technically redundant due to earlier filtering
        if np.isnan(lat) or np.isnan(lon) or np.isnan(ssh_norm):
            continue

        # Convert SSH norm value to a color
        # Apply colormap
        color = color_map(ssh_norm)
        folium.CircleMarker(
            location=[lat, lon],
            radius=1,  # Adjust the size as needed
            fill=True,
            color=color,
            fill_color=color,
            fill_opacity=0.7
        ).add_to(map_obj)


# Path to the NetCDF file
nc_file_path = 'Data/SWOT_L2_LR_SSH_Unsmoothed_555_015_20230617T211922_20230617T220935_PIB0_01.nc'

# Initialize the folium map
map = folium.Map(location=[0, 0], zoom_start=3, tiles='CartoDB positron')

# Process data in chunks
chunk_size = 1000  # Define a reasonable chunk size

try:
    with Dataset(nc_file_path, "r") as netcdf_file:
        # Assuming the total number of lines is known or can be calculated
        start = 19500
        end = 20000

        print("Plotting Left SSH")
        for start_line in range(start, end, chunk_size):
            end_line = min(start_line + chunk_size, end)

            # Read chunk data
            ssh_chunk = netcdf_file.groups['left'].variables['ssh_karin_2'][start_line:end_line, :].flatten()
            lat_chunk = np.ma.filled(netcdf_file.groups['left'].variables['latitude'][start_line:end_line, :].flatten(),
                                     np.nan)
            lon_chunk = np.ma.filled(
                netcdf_file.groups['left'].variables['longitude'][start_line:end_line, :].flatten(), np.nan)

            # Process and add this chunk's data to the map
            process_ssh_chunk(ssh_chunk, lat_chunk, lon_chunk, map)

            print(f"Left: Processing chunk from line {start_line} to {end_line}")

            print("Plotting Right SSH")
            for start_line in range(start, end, chunk_size):
                end_line = min(start_line + chunk_size, end)

                # Read chunk data
                ssh_chunk = netcdf_file.groups['right'].variables['ssh_karin_2'][start_line:end_line, :].flatten()
                lat_chunk = np.ma.filled(
                    netcdf_file.groups['right'].variables['latitude'][start_line:end_line, :].flatten(),
                    np.nan)
                lon_chunk = np.ma.filled(
                    netcdf_file.groups['right'].variables['longitude'][start_line:end_line, :].flatten(), np.nan)

                # Process and add this chunk's data to the map
                process_ssh_chunk(ssh_chunk, lat_chunk, lon_chunk, map)

                print(f"Right: Processing chunk from line {start_line} to {end_line}")

except Exception as e:
    print("An error occurred:", e)



# Save the map
print("Saving...")
map.save('ssh_map.html')

