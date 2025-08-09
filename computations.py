######CLASS COMPUTATIONS############
import scipy.signal as signal
from scipy.signal import welch, csd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from ncFileReader import RightSideData, LeftSideData

class Computations:

    #Bunch of constants
    g = 9.81 #Gravity [m/s^2]
    rho = 1000  # Water density [kg/m^3]
    #h = {'Shallow': 4000, 'Middle': 5000, 'Deep': 6000} #Ocean Depth, in m
    h = 5000

    # Min and Max Time Period of IG waves
    T_min = 60 #[sec]
    T_max = 200 #[sec]

    @staticmethod
    def dispersion_relation_f_from_k(k_prime):
        g = Computations.g
        h = Computations.h
        k = k_prime * 2 * np.pi
        """
        Calculate the frequency (f) from the wavelength (L) for given water depth (h).

        Parameters:
        L (float): Wavelength (m)
        h (float): Water depth (m)

        Returns:
        float: Frequency (Hz)
        """
        omega_squared = g * k * np.tanh(k * h)
        omega = np.sqrt(omega_squared)
        f = omega / (2 * np.pi)
        return f

    @staticmethod
    def dispersion_relation_L_from_f(f):
        g = Computations.g
        h = Computations.h
        """
        Calculate the wavelength (L) from the frequency (f) for given water depth (h).

        Parameters:
        f (float): Frequency (Hz)
        h (float): Water depth (m)

        Returns:
        float: Wavelength (m)
        """
        omega = 2 * np.pi * f
        k_initial = omega ** 2 / g

        # Iteratively solve for k using Newton-Raphson method
        def dispersion_relation(k):
            return g * k * np.tanh(k * h) - omega ** 2

        def dispersion_relation_prime(k):
            return g * np.tanh(k * h) + g * k * h * (1 / np.cosh(k * h)) ** 2

        k = k_initial
        tolerance = 1e-6
        max_iterations = 1000

        for _ in range(max_iterations):
            f_k = dispersion_relation(k)
            f_k_prime = dispersion_relation_prime(k)
            k_new = k - f_k / f_k_prime
            if np.abs(k_new - k) < tolerance:
                k = k_new
                break
            k = k_new

        L = 2 * np.pi / k
        return L

    @staticmethod
    def fill_in_nan(ssh_array):
        # fill in Nan values with ssh mean
        mean_ssh = np.ma.mean(ssh_array)
        new_ssh = np.ma.filled(ssh_array, mean_ssh)

        # normalize signal by mean
        new_ssh -= mean_ssh
        return new_ssh

    @staticmethod
    def sample_rate(ssh, theta):
        # calculating sample rate L [m]
        a = 60000
        if theta <= 45:
            rad_angle = theta / 180 * np.pi
            c = a / np.cos(rad_angle)

        elif  90 >= theta > 45:
            rad_angle = (90 - theta) / 180 * np.pi
            c = a / np.cos(rad_angle)

        elif  135 >= theta > 90:
            rad_angle = (90 - theta) / 180 * np.pi
            c = a / np.cos(rad_angle)

        else:
            rad_angle = (180 - theta) / 180 * np.pi
            c = a / np.cos(rad_angle)

        n = len(ssh)
        L = c / n
        return L


    @staticmethod
    # function plots the SSH values [m] for the x-array of angle theta
    # theta - angle in degrees
    def plot_x_vs_ssh(signal, theta, filter_type, difference):
        #calculate sample rate [m]
        L = Computations.sample_rate(signal, theta)

        plt.figure(figsize=(14, 10))
        signal_len = len(signal)
        start, end = 0, signal_len * L
        x = np.linspace(start, end, signal_len)

        plt.plot(x, signal)
        plt.scatter(x, signal, c='r', s=6)
        plt.gca().yaxis.set_major_locator(MaxNLocator(nbins='auto'))
        plt.xlabel("Distance [m]")
        plt.ylabel("SSH")

        if difference is False:
            plt.title(f"SSH of angle {theta} with {filter_type} Filter")
        if difference is True:
            plt.title(f'Difference between SSH of BP Filter and Downsampling (angle: {theta})')
        plt.show()


    @staticmethod
    def fft(ssh, theta):
        #calculate sample rate [metres]
        L = Computations.sample_rate(ssh, theta)

        # Calculate FFT
        fft_values = np.fft.fft(ssh)
        n = len(ssh)
        k_frequency = np.fft.fftfreq(n, d=L)    #K frq is (1/L)
        magnitude = np.abs(fft_values) / n #WE SHOULD NORMILIZE THIS BY n IN ORDER TO GET THE PYSICAL AMPLITUDES!
        phase = np.angle(fft_values)

        # Taking only positive frequencies
        positive_freq_indices = k_frequency > 0
        k_frequencies = k_frequency[positive_freq_indices]
        magnitudes = magnitude[positive_freq_indices]
        phases = phase[positive_freq_indices]

        return k_frequencies, magnitudes, phases

    @staticmethod
    # x - signal array
    # cutoff - cutoff freq (maximum: fs/2)
    def fft_with_lowpass(x, cutoff_freq, theta):
        #calculate sample rate [metres]
        L = Computations.sample_rate(x, theta)

        # applying low-pass butterworth filter
        fs = 1/L  # Sampling frequency
        nyq_freq = fs/2
        cutoff_rel = cutoff_freq/nyq_freq
        sos = signal.butter(3, cutoff_rel, 'low',  False, 'sos')
        filtered_signal = signal.sosfilt(sos, x)

        # Calculate FFT
        frequencies, magnitudes, phases = Computations.fft(filtered_signal, theta)
        return frequencies, magnitudes, phases, filtered_signal


    @staticmethod
    def plot_k_spectrum(k, mags, phases, theta, filter_type, difference):
        f = Computations.dispersion_relation_f_from_k(k)
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

        # Amplitude vs Wavenumber plot
        ax1.plot(k, mags, label='Magnitude')
        ax1.scatter(k, mags, c='r', s=6)
        ax1.set_xlabel("Wavenumber (1/m)")
        ax1.set_ylabel("SSH [m]")
        ax1.grid(True)

        ax1_tw = ax1.twiny()
        ax1_tw.set_xlim(ax1.get_xlim())
        ax1_tw.set_xticks(k[::len(k)//10])  # Select fewer ticks for clarity
        ax1_tw.set_xticklabels([f"{freq:.2e}" for freq in f[::len(f)//10]])
        ax1_tw.set_xlabel("Frequency (Hz)")

        # Phase vs Wavenumber plot
        mod_phase = np.mod(phases, 2 * np.pi)
        limit_phase = np.where(mod_phase > np.pi, mod_phase - 2 * np.pi, mod_phase)
        ax2.plot(k, limit_phase, label='Phase')
        ax2.scatter(k, limit_phase, c='r', s=6)
        ax2.set_xlabel("Wavenumber (1/m)")
        ax2.set_ylabel("Phase (rad)")
        ax2.grid(True)

        ax2_tw = ax2.twiny()
        ax2_tw.set_xlim(ax2.get_xlim())
        ax2_tw.set_xticks(k[::len(k)//10])  # Select fewer ticks for clarity
        ax2_tw.set_xticklabels([f"{freq:.2e}" for freq in f[::len(f)//10]])
        ax2_tw.set_xlabel("Frequency (Hz)")

        # Highlight the area between T=60sec and T=200sec
        f_min = 1 / Computations.T_min
        f_max = 1 / Computations.T_max

        # Find the corresponding k values for f_min and f_max
        k_min = 1/Computations.dispersion_relation_L_from_f(f_min)
        print(Computations.dispersion_relation_L_from_f(f_min))
        k_max = 1/Computations.dispersion_relation_L_from_f(f_max)
        print(Computations.dispersion_relation_L_from_f(f_max))

        ax1.axvspan(k_min, k_max, color='yellow', alpha=0.3)
        ax2.axvspan(k_min, k_max, color='yellow', alpha=0.3)

        if difference is False:
            plt.suptitle(f"Frequency Analysis of SSH: angle {theta} with {filter_type} filter")
        if difference is True:
            plt.suptitle(f'Frequency Domain Difference before and after {filter_type} filter')
        plt.tight_layout()
        plt.show()

    @staticmethod
    def bandpass(filtered_x, theta):
        #calculate sample rate [m]
        L = Computations.sample_rate(filtered_x, theta)

        # applying butterworth bandpass filter
        fs = 1 / L  # Sampling frequency
        nyq_freq = fs / 2
        f_low = 0.0000246522/nyq_freq
        f_high = 0.00017791/nyq_freq
        sos = signal.butter(3, [f_low, f_high], 'bandpass', False, 'sos')
        filtered_signal = signal.sosfilt(sos, filtered_x)

        # calculate FFT
        frequencies, magnitudes, phases = Computations.fft(filtered_signal, theta)

        return frequencies, magnitudes, phases, filtered_signal

    @staticmethod
    # function downsamples array by a factor "factor"
    def downsampling(x_array, factor):
        downsampled_x = x_array[::factor]

        return downsampled_x

    @staticmethod
    def find_max_freq(freq, amps):
        max_amps_ind, _ = signal.find_peaks(amps)
        max_freq = freq[max_amps_ind]
        ig_indexes = np.where((max_freq >= 0.0000246522) & (max_freq <= 0.00017791))
        max_freq = max_freq[ig_indexes]
        max_amps = amps[ig_indexes]

        return max_freq, max_amps

    @staticmethod
    def cross_correlation(ssh1, ssh2):
        if len(ssh1) != len(ssh2):
            raise ValueError("Input arrays must have the same length")

        mean_ssh1 = np.mean(ssh1)
        mean_ssh2 = np.mean(ssh2)
        std_ssh1 = np.std(ssh1)
        std_ssh2 = np.std(ssh2)

        if std_ssh1 == 0 or std_ssh2 == 0:
            return None  # Return None if standard deviation is zero

        correlation = np.correlate(ssh1 - mean_ssh1, ssh2 - mean_ssh2, mode='full')
        n = len(ssh1)
        return correlation / (n * std_ssh1 * std_ssh2)

    @staticmethod
    def plot_cross_correlation(ssh1, ssh2, theta):
        """
        Plot the cross-correlation between two SSH arrays with lag in meters.

        Parameters:
        ssh1 (array): First SSH array
        ssh2 (array): Second SSH array
        theta (float): Angle in degrees for sample rate calculation
        """
        L = Computations.sample_rate(ssh1, theta)
        correlation = Computations.cross_correlation(ssh1, ssh2)
        lags = np.arange(-len(ssh1) + 1, len(ssh1))
        lags_in_meters = lags * L

        plt.figure(figsize=(12, 6))
        plt.plot(lags_in_meters, correlation)
        plt.xlabel('Lag (meters)')
        plt.ylabel('Correlation')
        plt.title('Cross-Correlation between SSH Arrays')
        plt.grid(True)
        plt.show()

    @staticmethod
    def plot_multiple_cross_correlations(nc_file_path, line, theta_step=5):
        theta_values = np.arange(125, 126, theta_step)
        plt.figure(figsize=(12, 8))

        for theta in theta_values:
            print(theta)
            rightdata = RightSideData(nc_file_path)
            leftdata = LeftSideData(nc_file_path)

            rightSSH = rightdata.get_SSH_array_with_angle(line, theta)
            leftSSH = leftdata.get_SSH_array_with_angle(line, theta)

            _, _, _, rightSSH = Computations.bandpass(rightSSH, theta)
            _, _, _, leftSSH = Computations.bandpass(leftSSH, theta)

            correlation = Computations.cross_correlation(rightSSH, leftSSH)
            if correlation is None:
                print(f"Skipping theta = {theta}°: standard deviation of one of the arrays is zero")
                continue

            L = Computations.sample_rate(rightSSH, theta)
            lags = np.arange(-len(rightSSH) + 1, len(rightSSH))
            lags_in_meters = lags * L

            plt.plot(lags_in_meters, correlation, label=f'Theta = {theta}°')

        plt.xlabel('Lag (meters)')
        plt.ylabel('Correlation')
        plt.title('Cross-Correlation between SSH Arrays for Different Thetas')
        plt.legend(loc='upper right')
        plt.grid(True)
        plt.show()

    @staticmethod
    def cross_spectral_density(ssh1, ssh2, fs):
        """
        Compute the cross-spectral density between two SSH arrays.

        Parameters:
        ssh1 (array): First SSH array
        ssh2 (array): Second SSH array
        fs (float): Sampling frequency

        Returns:
        f (array): Frequency array
        Pxy (array): Cross-spectral density
        """
        f, Pxy = csd(ssh1, ssh2, fs=fs, nperseg=256)
        return f, Pxy

    @staticmethod
    def plot_cross_spectral_density(ssh1, ssh2, theta):
        """
        Plot the cross-spectral density between two SSH arrays.

        Parameters:
        ssh1 (array): First SSH array
        ssh2 (array): Second SSH array
        theta (float): Angle in degrees for sample rate calculation
        """
        L = Computations.sample_rate(ssh1, theta)
        fs = 1 / L  # Sampling frequency

        f, Pxy = Computations.cross_spectral_density(ssh1, ssh2, fs)

        plt.figure(figsize=(12, 6))

        # Plot magnitude of cross-spectral density
        plt.subplot(2, 1, 1)
        plt.semilogy(f, np.abs(Pxy))
        plt.xlabel('Frequency [Hz]')
        plt.ylabel('CSD Magnitude')
        plt.title('Cross-Spectral Density Magnitude')
        plt.grid(True)

        # Plot phase of cross-spectral density
        plt.subplot(2, 1, 2)
        plt.plot(f, np.angle(Pxy))
        plt.xlabel('Frequency [Hz]')
        plt.ylabel('CSD Phase [radians]')
        plt.title('Cross-Spectral Density Phase')
        plt.grid(True)

        plt.tight_layout()
        plt.show()

    @staticmethod
    def compute_energy_spectrum(ssh, theta):
        k_frq, magnitudes, phases = Computations.fft(ssh, theta)
        rho = Computations.rho
        g = Computations.g
        energy_spectrum = 0.5 * rho * g * magnitudes ** 2
        return k_frq, energy_spectrum, phases