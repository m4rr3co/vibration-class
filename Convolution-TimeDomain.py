# Vibrations Course
# Mechanical Engineering Post Graduation Program
# UNESP - Ilha Solteira, SP - Brazil
# ####
# Work by Lucas Veronez Goulart Ferreira
# ####

import numpy as np
import matplotlib.pyplot as plot
from scipy import signal

# Time domain convolution: f(t)*h(t) = x(t)
# 1st case: f(t) = delayed impulse

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

sampling_frequency = 500
t_max = 10
time_range = t = np.arange(0,t_max,1/sampling_frequency)

h = (np.exp(-damping_ratio*wn*t)/wd*m)*np.sin(wd*t)
imp = signal.unit_impulse(len(time_range),int(len(time_range)*0.1))

conv = np.convolve(h,imp)

plot.plot(conv)
plot.show()
