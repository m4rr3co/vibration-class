# Vibrations Course
# Mechanical Engineering Post Graduation Program
# UNESP - Ilha Solteira, SP - Brazil
# ####
# Work by Lucas Veronez Goulart Ferreira
# ####

import numpy as np
import matplotlib.pyplot as plot
from scipy.signal import hilbert, chirp

# Forced Vibration of a Single-Degree of Freedom Mass-Spring-Dampener System
# Objective: Frequency Response Function
# Solution: H(wj) = X/F = 1 / (k-m*w^2+j*w*c)
# Parameters
mass = m = 1  # Kilograms [kg]
stiffness = 10e4  # Newton/meter [N/m]
wn = np.sqrt(stiffness/mass)  # Natural frequency
critical_damping = cc = 2*m*wn  # zeta = c/cc

# Damping Ratio
# For now, these are given, but can be calculated from
# response as will be shown further on.
damping_ratio = [0.1,0.01,0.001]  # No dimension, i.e.: zeta = c/sqrt(m.k)
c = np.array(damping_ratio)*cc  # Damping values
print('Damping ratios (Zeta): '+str(damping_ratio))
print('Damping values [N.s/m]: '+str(c))
pi = np.pi  # Pi = Everyone knows what this is
interval = t = np.arange(0,1.5,0.001)  # Time interval [0,0.001,...,10] seconds


colors = ['r','g','b']


# for i in range(0,3):
    # x = (np.exp(-damping_ratio[i]*wn*t)/wd[i]*m)*np.sin(wd[i]*t)
    # plot.plot(t,x,c[i])
    # plot.plot(t,x,c[i])
#plot.show()
