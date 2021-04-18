# Vibrations Course
# Mechanical Engineering Post Graduation Program
# UNESP - Ilha Solteira, SP - Brazil
# ####
# Work by Lucas Veronez Goulart Ferreira
# ####

import numpy as np
import matplotlib.pyplot as plot
from scipy.fft import ifft

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

# Sampling Parameters
frequency_range = 6*fn
sample_spacing = 0.01
number_of_samples = int(frequency_range/sample_spacing)

#  freq_bins: The X axis represents the frequency components of the sampled signals in Hz.
freq_bins = f = np.linspace(-frequency_range/2,frequency_range/2,number_of_samples)
frf = 1/(k-(m*((2 * np.pi * f)**2))+(1j * c * (2 * np.pi * f)))
ift = ifft(frf)/sample_spacing

plot.subplot(2,1,1)
plot.plot(freq_bins,np.abs(frf))

plot.subplot(2,1,2)
plot.plot(np.real(ift))
plot.show()
