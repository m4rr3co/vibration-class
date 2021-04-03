# Vibrations Course
# Mechanical Engineering Post Graduation Program
# UNESP - Ilha Solteira, SP - Brazil
# ####
# Work by Lucas Veronez Goulart Ferreira
# ####

import numpy as np
import matplotlib.pyplot as plot
from scipy.signal import hilbert


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


# Estimation of Damping from FRF
# FRF: H(wj) = X/F = 1 / (k-m*w^2+j*w*c)
# Parameters
mass = m = 1  # Kilograms [kg]
stiffness = k = 10e3  # Newton/meter [N/m]
wn = np.sqrt(k/m)  # Natural frequency

# Damping Ratio
# For now, these are given, but can be calculated from
# response as will be shown further on.
damping_ratio = [0.1,0.01,0.001]  # No dimension, i.e.: zeta = c/2*sqrt(m.k)
damping = c = np.asarray(damping_ratio)*2*np.sqrt(m*k)
pi = np.pi  # Pi = Everybody knows what this is
interval = w = np.arange(0,300,0.01)  # Time interval [0,0.001,...,100] rad/s

# Calculated Parameters
wd = [0,0,0]  # Defining array of values for Damped natural frequencies
for i in range(len(damping_ratio)):
    wd[i] = wn*np.sqrt(1-damping_ratio[i]**2)
print('Damped natural frequencies [rad/s]: '+str(wd))

colors = ['r','g','b']
cdr = [0,0,0]
j = 1j
for i in range(0,3):
    Hwj = 1 / ((k - (m * w**2)) + (w * c[i] * j))
    abs_of_Hwj = np.abs(Hwj)
    max_value = max(abs_of_Hwj)
    index_of_max = np.where(abs_of_Hwj == max_value)[0][0]
    first_part = abs_of_Hwj[:index_of_max]
    a = find_nearest(first_part,max_value/np.sqrt(2))
    second_part = abs_of_Hwj[index_of_max:]
    b = find_nearest(second_part, max_value / np.sqrt(2))
    index_of_3dblower1 = np.where(abs_of_Hwj == a)[0][0]
    index_of_3dblower2 = np.where(abs_of_Hwj == b)[0][0]
    cdr[i] = (w[index_of_3dblower2]-w[index_of_3dblower1])/(2*w[index_of_max])
print('Damping ratios (Zeta): '+str(damping_ratio))
print('Calculated damping ratios (Zeta): '+str(cdr))





