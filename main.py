import xarray as xr
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Path to the NetCDF file
nc_file_path = 'Data/SWOT_L2_LR_SSH_Unsmoothed_555_015_20230617T211922_20230617T220935_PIB0_01.nc'  # Update this path to your file's location


# Open the NetCDF file using xarray
ds = xr.open_dataset(nc_file_path)
print(ds)

# Specify the indices for the pixel you're interested in
line_index = 3  # For example, the first line
pixel_index = 0  # For example, the first pixel

ssh_karin_value = ds['ssh_karin'].isel(num_lines=line_index, num_pixels=pixel_index).values
latitude_value = ds['latitude'].isel(num_lines=line_index, num_pixels=pixel_index).values
longitude_value = ds['longitude'].isel(num_lines=line_index, num_pixels=pixel_index).values
time_value = ds['time'].isel(num_lines=line_index).values  # Time is typically associated with the line rather than individual pixels
# Convert the time_value to a more readable format if necessary
# Assuming time_value is in a standard format recognized by pandas/numpy, such as numpy.datetime64 or a similar
time_value_readable = pd.to_datetime(str(time_value))
time_value_readable = time_value_readable.strftime('%Y-%m-%d %H:%M:%S.%f')

# Display the properties of the pixel including the time tag
print(f"Time: {time_value_readable}, Latitude: {latitude_value}, Longitude: {longitude_value}, SSH (ssh_karin): {ssh_karin_value}")

# Don't forget to close the dataset after you're done
ds.close()