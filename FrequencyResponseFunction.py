# Vibrations Course
# Mechanical Engineering Post Graduation Program
# UNESP - Ilha Solteira, SP - Brazil
# ####
# Work by Lucas Veronez Goulart Ferreira
# ####

import numpy as np
import matplotlib.pyplot as plot

# Forced Vibration of a Single-Degree of Freedom Mass-Spring-Dampener System
# Objective: Frequency Response Function
# Solution: H(wj) = X/F = 1 / (k-m*w^2+j*w*c)
# Parameters
mass = m = 1  # Kilograms [kg]
stiffness = k = 10e4  # Newton/meter [N/m]
wn = np.sqrt(stiffness/mass)  # Natural frequency
critical_damping = cc = 2*m*wn  # zeta = c/cc

# Damping Ratio
# For now, these are given, but can be calculated from
# response as will be shown further on.
damping_ratio = [0.999,0.5,0.3]  # No dimension, i.e.: zeta = c/sqrt(m.k)
c = np.array(damping_ratio)*cc  # Damping values
print('Damping ratios (Zeta): '+str(damping_ratio))
print('Damping values [N.s/m]: '+str(c))
pi = np.pi  # Pi = Everyone knows what this is
interval = w = np.arange(0,50,0.01)  # Frequency interval [0,0.1,...,100] rad/s

colors = ['r','g','b']
j = 1j

for i in range(0,3):
    x = 1 / (k-(m*w**2)+(w*c[i]*j))
    imaginary_part = x.imag
    real_part = x.real
    phase = np.angle(x)
    # plot.plot(w,np.abs(x),colors[i])  # Absolute value
    # plot.plot(w,phase,colors[i])  # Phase plot
    # plot.plot(w,real_part,colors[i])  # Real part plot
    # plot.plot(w,imaginary_part,colors[i])  # Imaginary part plot
    # plot.plot(real_part,imaginary_part,colors[i])  # Nyquist plot
plot.show()
