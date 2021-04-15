# Vibrations Course
# Mechanical Engineering Post Graduation Program
# UNESP - Ilha Solteira, SP - Brazil
# ####
# Work by Lucas Veronez Goulart Ferreira
# ####

import numpy as np
import matplotlib.pyplot as plot
from scipy.fft import fft, fftfreq

# Sampling Parameters
sample_frequency = 120  # Hz
window = 12.58204  # seconds
start_time = 5.12837  # seconds
end_time = start_time+window  # seconds
number_of_samples = N = int(sample_frequency * window)  # seconds
time_interval = t = np.linspace(start_time,end_time,number_of_samples)

# Signal: Impulse Response Function
mass = m = 1  # [Kg]
stiffness = k = 10e3  # [N/m]
wn = np.sqrt(k/m)  # [rad/s] - Natural Frequency
fn = wn / (2 * np.pi)  # [Hz] - Natural Frequency
damping_ratio = 0.01  # Zeta
critical_damping = cc = 2*np.sqrt(k*m)  # [Ns/m]
c = damping_ratio * critical_damping  # [Ns/m] - Damping
wd = wn * np.sqrt(1-damping_ratio**2)  # [rad/s] - Damped Natural Frequency
fd = wd / (2 * np.pi)  # [Hz] - Damped Natural Frequency
print('Calculated Damped Natural Frequency = '+str(fd)+' Hz.')

# Impulse response function as a function of time in seconds.

#  freq_bins: The X axis represents the frequency components of the sampled signals in Hz.
freq_bins = f = np.linspace(0,sample_frequency,number_of_samples)

# Frequencies used here are calculated in Hz.
samples = y = (np.exp(-damping_ratio*(2 * np.pi * fn)*t)/(2 * np.pi * fd)*m)*np.sin(2 * np.pi * fd * t)
ft = fft(y)  # Fast Fourier Transform of the sampled signal

frf = 1/(k-(m*((2 * np.pi * f)**2))+(1j * c * (2 * np.pi * f)))

plot.ylabel('H')
plot.xlabel('[Hz]')
plot.plot(freq_bins,np.abs(ft),freq_bins,np.abs(frf))  # fft versus frequency bins chart
plot.show()

max_index = np.where(np.abs(ft) == max(np.abs(ft)))
print('Damped Natural Frequency from FFT: '+str(freq_bins[max_index])+' Hz.')
