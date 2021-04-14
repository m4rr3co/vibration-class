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
j = 0+1j  # Imaginary unit
mass = m = 1  # Kilograms [kg]
stiffness = k = 10e3  # Newton/meter [N/m]
damping_ratio = [0.1,0.01,0.001]
natural_frequency = wn = np.sqrt(stiffness/mass)
fn = wn / (2 * np.pi)
critical_damping = cc = 2*m*wn
c = cc * damping_ratio[0]
wd = wn*np.sqrt(1-(damping_ratio[0]**2))
fd = wd / (2 * np.pi)

angular_velocity = w = np.arange(0,300,1)
H_wj = 1 / (k-(m*w**2)+(w*c*j))

sample_spacing = dt = 0.001
window = 100
time_interval = t = np.arange(0,window,dt)
x = (np.exp(-damping_ratio[0]*(2 * np.pi * fn)*t)/(2 * np.pi * fd)*m)*np.sin((2 * np.pi * fd)*t)

H_dwj = fft(x)
freq_bins = np.linspace(0,window * 10,len(H_dwj))
plot.plot(fftfreq(len(H_dwj),dt),np.abs(H_dwj))
plot.show()
