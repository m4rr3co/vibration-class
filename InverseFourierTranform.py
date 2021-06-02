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
stiffness = k = 1e4  # [N/m]
wn = np.sqrt(k/m)  # [rad/s] - Natural Frequency
fn = wn / (2 * np.pi)  # [Hz] - Natural Frequency
damping_ratio = 0.1  # Zeta
critical_damping = cc = 2*np.sqrt(k*m)  # [Ns/m]
c = damping_ratio * critical_damping  # [Ns/m] - Damping
wd = wn * np.sqrt(1-damping_ratio**2)  # [rad/s] - Damped Natural Frequency
fd = wd / (2 * np.pi)  # [Hz] - Damped Natural Frequency
print('Calculated Damped Natural Frequency = '+str(fd)+' Hz.')

# Sampling Parameters
delta = 0.001  # seconds
sample_frequency = 1/delta  # Hz
tmax = 133.73  # seconds
time_interval = t = np.arange(0,tmax,delta)

#  freq_bins: The X axis represents the frequency components of the sampled signals in Hz.
freq_bins = f = np.linspace(0,sample_frequency,len(t))
frf = 1/(k-(m*((2 * np.pi * f)**2))+(1j * c * (2 * np.pi * f)))
double_sided_frf = ds_frf = frf+np.conj(np.flip(frf))

x_t = ifft(double_sided_frf)/delta

x_t_true = (np.exp(-damping_ratio*(2 * np.pi * fn)*t)/(2 * np.pi * fd)*m)*np.sin(2 * np.pi * fd * t)

plot.plot(t,np.real(x_t))
plot.plot(t,x_t_true)
plot.show()
