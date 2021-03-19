import numpy as np
import matplotlib.pyplot as plot

# Vibrations Course
# Mechanical Engineering Post Graduation Program
# UNESP - Ilha Solteira, SP - Brazil
# ####
# Work by Lucas Veronez Goulart Ferreira
# ####

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
k = stiffness = 1  # Newton/meter [N/m]
damping_ratio = 0.1

# Time domain analysis
critical_damping = 2*np.sqrt(stiffness*mass)
c = actual_damping = damping_ratio*critical_damping
Xo = f/np.sqrt(np.float_power(k-m*np.float_power(w,2),2)+np.float_power(w*c,2))
print('Amplitude of response = '+str(Xo)+' N')
phi = np.arctan(w*c/k-np.float_power(w,2)*m)
print('Phase angle = '+str(phi)+' rad')
x = Xo*np.sin(w*interval+phi)
# Plots
# plot.title('Red Curve - f(t) [N] \n Green Curve - response [mm]')
# plot.plot(interval,f_of_t,'r')
# plot.plot(interval,1000*x,'g')
# plot.show()

# Frequency domain analysis
# Euler's Equation: e^j.teta = cos(teta)+j.sin(teta)
# Euler's Equation: e^-j.teta = cos(teta)-j.sin(teta)
# Equation of motion: mx'' + cx' + kx = F . e^j.w.t
# Solution: x = X . e^j.w.t


# Forced Vibration of a Single-Degree of Freedom Mass-Spring-Dampener System
# Input: Impulse
# Parameters
mass = 1  # Kilograms [kg]
stiffness = 10e4  # Newton/meter [N/m]

# Damping Ratio
# For now, these are given, but can be calculated
# from response as will be shown further on.
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




