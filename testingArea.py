import math

from netCDF4 import Dataset
from ncFileReader import RightSideData, LeftSideData
from computations import Computations
import numpy as np
from dataVisualizer import DataVisualizer
'''
nc_file_path = 'Data/SWOT_L2_LR_SSH_Unsmoothed_555_015_20230617T211922_20230617T220935_PIB0_01.nc'
theta = 0
line = 45000

rightdata = RightSideData(nc_file_path)
leftdata = LeftSideData(nc_file_path)


rightSSH = rightdata.get_SSH_array_with_angle(line, theta)
leftSSH = leftdata.get_SSH_array_with_angle(line, theta)


_, _, _ ,rightSSH = Computations.bandpass(rightSSH, theta)
_, _, _ ,leftSSH = Computations.bandpass(leftSSH, theta)


Computations.plot_x_vs_ssh(rightSSH, theta, 'BP', None)
Computations.plot_x_vs_ssh(leftSSH, theta, 'BP', None)

Computations.plot_cross_correlation(rightSSH, leftSSH, theta)
Computations.plot_cross_spectral_density(rightSSH, leftSSH, theta)



# Example usage
h = 3500  # Water depth in meters

# Calculate frequency from wavelength
L = 10000  # Wavelength in meters
f = Computations.dispersion_relation_f_from_L(L, h)
print(f"Frequency (f) for L={L} and h={h}: {f:.6f} Hz, the time period is: {1/f}s")


# Calculate wavelength from frequency
f = 1/300  # Frequency in Hz
L = Computations.dispersion_relation_L_from_f(f, h)
print(f"Max Wavelength (L) for f={f} or time period: {1/f}s, and h={h}: {L:.6f} m")

# Calculate wavelength from frequency
f = 1/30  # Frequency in Hz
L = Computations.dispersion_relation_L_from_f(f, h)
print(f"Min Wavelength (L) for f={f} or time period: {1/f}s, and h={h}: {L:.6f} m")




'''

nc_file_path = 'Data/SWOT_L2_LR_SSH_Unsmoothed_001_377_20230803T155753_20230803T164920_PGC0_01.nc'

app = DataVisualizer(nc_file_path)
app.mainloop()

'''

k = 0.00017
print(1/Computations.dispersion_relation_f_from_k(k))


nc_file_path = 'Data/SWOT_L2_LR_SSH_Unsmoothed_555_015_20230617T211922_20230617T220935_PIB0_01.nc'
line = 47000
Computations.plot_multiple_cross_correlations(nc_file_path, line)
'''