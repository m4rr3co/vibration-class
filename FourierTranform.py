# Vibrations Course
# Mechanical Engineering Post Graduation Program
# UNESP - Ilha Solteira, SP - Brazil
# ####
# Work by Lucas Veronez Goulart Ferreira
# ####

import numpy as np
import matplotlib.pyplot as plot
from scipy.fft import fft

# Sampling Parameters
delta = 0.01  # seconds
sample_frequency = 1/delta  # Hz
tmax = 11.73  # seconds
time_interval = t = np.arange(0,tmax,delta)

# Signal: Impulse Response Function
mass = m = 1  # [Kg]
stiffness = k = 1e4  # [N/m]
wn = np.sqrt(k/m)  # [rad/s] - Natural Frequency
fn = wn / (2 * np.pi)  # [Hz] - Natural Frequency
damping_ratio = 0.1  # Zeta
critical_damping = cc = 2*np.sqrt(k*m)  # [Ns/m]
c = damping_ratio * critical_damping  # [Ns/m] - Damping
wd = wn * np.sqrt(1-(damping_ratio**2))  # [rad/s] - Damped Natural Frequency
fd = wd / (2 * np.pi)  # [Hz] - Damped Natural Frequency
print('Calculated Damped Natural Frequency = '+str(fd)+' Hz.')

# Frequencies used here are calculated in Hz.
f1 = np.linspace(0,sample_frequency,len(t))
samples = y = (np.exp(-damping_ratio*(2 * np.pi * fn)*t)/(2 * np.pi * fd)*m)*np.sin(2 * np.pi * fd * t)
ft = fft(y)*delta  # Fast Fourier Transform of the sampled signal

#  freq_bins: The X axis represents the frequency components of the sampled signal in Hz.
freq_bins = f = np.linspace(0,sample_frequency,len(t))
frf = 1/(k-(m*((2 * np.pi * f)**2))+(1j * c * (2 * np.pi * f)))

plot.ylabel('X/F = H')
plot.xlabel('Frequency [Hz]')
plot.plot(f,np.abs(frf),label='Analytical FRF')  # frf function versus frequency bins
plot.plot(f,np.abs(ft),label='FFT of h(t)')  # fft versus frequency bins chart
plot.legend()
plot.grid()
plot.show()

max_index = np.where(np.abs(ft) == max(np.abs(ft)))
print('Damped Natural Frequency from FFT: '+str(freq_bins[max_index])+' Hz.')
