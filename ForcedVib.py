# Vibrations Course
# Mechanical Engineering Post Graduation Program
# UNESP - Ilha Solteira, SP - Brazil
# ####
# Work by Lucas Veronez Goulart Ferreira
# ####

import numpy as np
import matplotlib.pyplot as plot

# Forced Vibration of a Single-Degree of Freedom
# Input: Harmonic Force f(t) = F.sin(w.t)
# Equation of motion: mx'' + cx' + kx = F.sin(w.t)
# Solution: x = Xo . sin(w.t + phi)

# Example: F = 10N, w = 10pi rad/s, m = 1kfg, Damping ratio = 0.1, k = 10 N/m
interval = np.arange(0,2,0.001)  # Time interval [0,0.001,...,10] seconds
f = force = 10  # Newtons
w = omega = 10*np.pi  # rad/s
f_of_t = f*np.sin(w*interval)
m = mass = 1  # Kilograms [kg]
k = stiffness = 10e4  # Newton/meter [N/m]
damping_ratio = 0.1

# Time domain analysis
critical_damping = 2*np.sqrt(stiffness*mass)
c = actual_damping = damping_ratio*critical_damping
Xo = f/np.sqrt(np.float_power(k-m*np.float_power(w,2),2)+np.float_power(w*c,2))
print('Amplitude of response = '+str(Xo)+' N')
phi = np.arctan(w*c/k-np.float_power(w,2)*m)
print('Phase angle = '+str(phi)+' rad')
x = Xo*np.sin(w*interval+phi)
deltast = f/k  # Deflection under static force f
# Plots
# plot.title('Red Curve - f(t) [N] \n Green Curve - response [mm]')
# plot.plot(interval,f_of_t,'r')
# plot.plot(interval,1000*x,'g')
# plot.show()

# zeta = damping_ratios = np.concatenate((np.arange(0.1,1,0.1),np.arange(1,5,0.5)))
# for damping_ratio in damping_ratios:
#    r = frequency_ratios = np.arange(0,5,0.01)
#    amplitude_over_deltast = 1/np.sqrt(np.float_power((1-np.float_power(r,2)),2)+np.float_power(2*damping_ratio*r,2))
#    plot.plot(frequency_ratios,amplitude_over_deltast)
# plot.show()

# Frequency domain analysis
# Euler's Equation: e^j.teta = cos(teta)+j.sin(teta)
# Euler's Equation: e^-j.teta = cos(teta)-j.sin(teta)
# Equation of motion: mx'' + cx' + kx = F . e^j.w.t
# Solution: x = X . e^j.w.t
