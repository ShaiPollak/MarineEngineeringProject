import numpy as np
from computations import Computations
from ncFileReader import RightSideData, LeftSideData
import matplotlib.pyplot as plt

# Path to the NetCDF file
nc_file_path = 'Data/SWOT_L2_LR_SSH_Unsmoothed_001_377_20230803T155753_20230803T164920_PGC0_01.nc'
rightdata = RightSideData(nc_file_path)

def find_ks():
    thetas = [*range(0, 180, 5)]
    all_frequencies = []
    all_thetas = []
    all_lines = []


    for idx, line in enumerate(range(45000, 50000, 1000)):
        frequencies = []
        t = []
        for theta in thetas:
            new_ssh = rightdata.get_SSH_array_with_angle(line, theta)
            #new_ssh = rightdata.get_SSH_line(52000)

            # Plot the SSH values of the 1d array
            #Computations.plot_x_vs_ssh(new_ssh, theta, 'No', difference=False)

            # Save the Frequency Domain values for the SSH array with no filter
            ssh_freqs, ssh_amps, ssh_phase = Computations.fft(new_ssh, theta)

            # Save the Frequency Domain values for the SSH array with lowpass filter
            freqs_lp, amps_lp, phase_lp, filtered_sig = Computations.fft_with_lowpass(new_ssh, 0.00017791, theta)

            # Plot the SSH values after the LP filter
            #Computations.plot_x_vs_ssh(filtered_sig,  theta, 'LP', difference=False)

            # Plot all the K spectrum Plots
            #Computations.plot_k_spectrum(ssh_freqs, ssh_amps, ssh_phase, theta, "No", difference=False)
            #Computations.plot_k_spectrum(freqs_lp, amps_lp, phase_lp, theta, 'LP', difference=False)

            diff_amps = amps_lp - ssh_amps
            diff_phases = phase_lp - ssh_phase
            #Computations.plot_k_spectrum(freqs_lp, diff_amps, diff_phases, theta, 'LP', difference=True)

            # Applying Bandpass Filter to Signal
            freqs_bp, mags_bp, phase_bp, filtered_BP_sig = Computations.bandpass(filtered_sig, theta)

            # Plot SSH for BP
            #Computations.plot_x_vs_ssh(filtered_BP_sig, theta, 'LP & BP', difference=False)

            # Plot K spectrum plots for BP
            #Computations.plot_k_spectrum(freqs_bp, mags_bp, phase_bp, 0, 'LP & BP', difference=False)

            # Find K amplitude peaks in the frequency domain of the BP
            max_freqs, max_amps = Computations.find_max_freq(freqs_bp, mags_bp)

            frequencies.extend(max_freqs)
            n = len(max_amps)
            for i in range(n):
                t.append(theta)

        all_frequencies.extend(frequencies)
        all_thetas.extend(t)
        all_lines.extend([idx] * len(frequencies))

    plt.figure(figsize=(10, 6))
    scatter = plt.scatter(all_thetas, all_frequencies, c=all_lines, cmap='jet', alpha=0.6)
    plt.title('Find Ks')
    plt.xlabel('Theta')
    plt.ylabel('K')
    plt.colorbar(scatter, ticks=range(len(range(45000, 50000, 1000))), label='Lines')
    plt.grid(True)
    plt.show()


find_ks()