import numpy as np
from computations import Computations
from ncFileReader import RightSideData, LeftSideData

# Path to the NetCDF file
nc_file_path = 'Data/SWOT_L2_LR_SSH_Unsmoothed_001_377_20230803T155753_20230803T164920_PGC0_01.nc'
theta = 125

rightdata = RightSideData(nc_file_path)
new_ssh = rightdata.get_SSH_array_with_angle(45000, theta)
#new_ssh = rightdata.get_SSH_line(52000)

# Plot the SSH values of the 1d array
Computations.plot_x_vs_ssh(new_ssh, theta, 'No', difference=False)

# Save the Frequency Domain values for the SSH array with no filter
ssh_freqs, ssh_amps, ssh_phase = Computations.compute_energy_spectrum(new_ssh, theta)

# Save the Frequency Domain values for the SSH array with lowpass filter
freqs_lp, amps_lp, phase_lp, filtered_sig = Computations.fft_with_lowpass(new_ssh, 0.00017791, theta)

# Plot the SSH values after the LP filter
Computations.plot_x_vs_ssh(filtered_sig,  theta, 'LP', difference=False)

# Plot all the K spectrum Plots
Computations.plot_k_spectrum(ssh_freqs, ssh_amps, ssh_phase, theta, "No", difference=False)
Computations.plot_k_spectrum(freqs_lp, amps_lp, phase_lp, theta, 'LP', difference=False)

diff_amps = amps_lp - ssh_amps
diff_phases = phase_lp - ssh_phase
Computations.plot_k_spectrum(freqs_lp, diff_amps, diff_phases, theta, 'LP', difference=True)

# Applying Bandpass Filter to Signal
freqs_bp, mags_bp, phase_bp, filtered_BP_sig = Computations.bandpass(filtered_sig, theta)


# Downsample SSH by a factor of 20
ds_ssh = Computations.downsampling(filtered_sig, 10)
ds_freqs, ds_amps, ds_phases = Computations.fft(ds_ssh, theta)

# Plot SSH for BP and DS
Computations.plot_x_vs_ssh(filtered_BP_sig, theta, 'LP & BP', difference=False)
Computations.plot_x_vs_ssh(ds_ssh, theta, "LP and Downsampling", difference=False)


# Plot K spectrum plots for BP and DS signals
Computations.plot_k_spectrum(freqs_bp, mags_bp, phase_bp, 0, 'LP & BP', difference=False)
Computations.plot_k_spectrum(ds_freqs, ds_amps, ds_phases, 0, "Downsampled, No", difference=False)

#Comparing BP filtered single with DS signal
freqs_bp_compare = freqs_bp[:len(ds_freqs)]
mags_bp_compare = mags_bp[:len(ds_amps)]
phase_bp_compare = phase_bp[:len(ds_phases)]


Computations.plot_k_spectrum(freqs_bp_compare, mags_bp_compare, phase_bp_compare, theta, 'BP for comparison', difference=False)
amps_diff = mags_bp_compare - ds_amps
phase_diff = phase_bp_compare - ds_phases

Computations.plot_k_spectrum(ds_freqs, amps_diff, phase_diff, theta, "BP and Downsampling", difference=True)

# Find K amplitude peaks in the frequency domain of the BP and DS signals
max_freqs, max_amps = Computations.find_max_freq(freqs_bp, mags_bp)
max_freqs_ds, max_amps_ds = Computations.find_max_freq(ds_freqs, ds_amps)

print(f'Bandpass filter max frequencies: {max_freqs} and amplitudes: {max_amps}')
print(f'Downsampled signal max frequencies: {max_freqs_ds} and amplitudes: {max_amps_ds}')


