# Vibrations Course
# Mechanical Engineering Post Graduation Program
# UNESP - Ilha Solteira, SP - Brazil
# ####
# Work by Lucas Veronez Goulart Ferreira
# ####

import numpy as np
import matplotlib.pyplot as plot
from scipy.fft import fft, fftfreq

# Free Vibration of a Single-Degree of Freedom Mass-Spring-Dampener System
# Input: Impulse
# Output: x(t) = (e^(zeta*wn*t)/m*wd)*sin(wd*t)
# Parameters
mass = m = 1  # Kilograms [kg]
stiffness = k = 10e3  # Newton/meter [N/m]
damping_ratio = [0.1,0.01,0.001]

# Imaginary unit
j = 0+1j

# Calculated Parameters
natural_frequency = wn = np.sqrt(stiffness/mass)
critical_damping = cc = 2*m*wn
wd = [0,0,0]  # Defining array of values for Damped natural frequencies
damping = c = [0,0,0]
damped_frequency = df = [0,0,0]
for i in range(len(damping_ratio)):
    wd[i] = wn*np.sqrt(1-np.float_power(damping_ratio[i],2))
    c[i] = damping_ratio[i]*cc
    df[i] = wd[i]/(2*np.pi)

print('Natural frequency: '+str(wn)+' rad/s.')
print('Damping ratios (Zeta): '+str(damping_ratio))
print('Damped natural frequencies [rad/s]: '+str(wd))
print('Damped natural frequencies [Hz]: '+str(df))

sample_spacing = 0.001
interval = t = np.arange(0,0.3,sample_spacing)  # Time interval [0,0.001,0.002,...] seconds
x = (np.exp(-damping_ratio[0]*wn*t)/(wd[0]*m)) * np.sin(wd[0]*t)
fourier = fft(x)
freqs = f = fftfreq(len(x),sample_spacing)
xf = (1/(k-(m*((2*np.pi*f)**2))+((2*np.pi*f)*c[0]*j)))
x1 = np.fft.ifft(xf)
print(len(x))
print(len(fourier))
print(len(freqs))

# plot.plot(t,x1)
# plot.plot(t,x)
# plot.plot(freqs,np.abs(xf))
plot.plot(freqs,np.abs(fourier))
plot.show()
