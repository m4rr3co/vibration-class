# Vibrations Course
# Mechanical Engineering Post Graduation Program
# UNESP - Ilha Solteira, SP - Brazil
# ####
# Work by Lucas Veronez Goulart Ferreira
# ####

import numpy as np
import matplotlib.pyplot as plot

# Forced Vibration of a Single-Degree of Freedom Mass-Spring-Dampener System
# Input: Impulse
# Parameters
mass = 1  # Kilograms [kg]
stiffness = 10e4  # Newton/meter [N/m]

# Damping Ratio
# For now, these are given, but can be calculated from
# response as will be shown further on.
damping_ratio = [0.1,0.01,0.001]  # No dimension, i.e.: zeta = c/sqrt(m.k)
pi = np.pi  # Pi = Everyone knows what this is
interval = np.arange(0,10,0.001)  # Time interval [0,0.001,...,10] seconds

# Calculated Parameters
wn = np.sqrt(stiffness/mass)  # Natural frequency
print('Natural frequency = '+str(wn)+' Hz')

wd = [0,0,0]  # Defining array of values for Damped natural frequencies
for i in range(len(damping_ratio)):
    wd[i] = wn*np.sqrt(1-np.float_power(damping_ratio[i],2))
print('Damping ratios (Zeta): '+str(damping_ratio))
print('Damped natural frequencies [Hz]: '+str(wd))




