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
stiffness = k = 1e4  # Newton/meter [N/m]
wn = np.sqrt(stiffness/mass)  # Natural frequency
critical_damping = cc = 2*m*wn  # zeta = c/cc

# Damping Ratio
# For now, these are given, but can be calculated from
# response as will be shown further on.
damping_ratio = [0.1,0.01,0.001]  # No dimension, i.e.: zeta = c/sqrt(m.k)
c = np.array(damping_ratio)*cc  # Damping values
frequency_interval = f = np.arange(0,30,0.01)  # Frequency interval [Hz]

colors = ['r','g','b']

for i in range(0,3):
    x = 1 / (k-(m*(2*np.pi*f)**2)+((2*np.pi*f)*c[i]*1j))
    imaginary_part = x.imag
    real_part = x.real
    phase = np.angle(x)
    plot.plot(f,np.abs(x),colors[i],label='\u03B6 = '+str(damping_ratio[i]))  # Absolute value
    plot.ylabel('X/F [m/N]')
    plot.xlabel('Frequency [Hz]')
    plot.yscale('log')
    # plot.plot(w,phase,colors[i])  # Phase plot
    # plot.plot(w,real_part,colors[i])  # Real part plot
    # plot.plot(w,imaginary_part,colors[i])  # Imaginary part plot
    # plot.plot(real_part,imaginary_part,colors[i])  # Nyquist plot
plot.axvline(x=(wn/(2*np.pi)),linestyle='--')
plot.legend()
plot.show()

for i in range(0,3):
    x = 1 / (k-(m*(2*np.pi*f)**2)+((2*np.pi*f)*c[i]*1j))
    imaginary_part = x.imag
    real_part = x.real
    phase = np.angle(x)
    # plot.plot(f,np.abs(x),colors[i],label='\u03B6 = '+str(damping_ratio[i]))  # Absolute value
    plot.subplot(2,2,1)
    plot.ylabel('Phase \u03C6 [rad]')
    plot.xlabel('Frequency [Hz]')
    plot.plot(f,phase,colors[i],label='\u03B6 = '+str(damping_ratio[i]))  # Phase plot

    plot.subplot(2, 2, 2)
    plot.ylabel('Real part [m/N]')
    plot.xlabel('Frequency [Hz]')
    plot.plot(f,real_part,colors[i],label='\u03B6 = '+str(damping_ratio[i]))  # Real part plot

    plot.subplot(2, 2, 3)
    plot.ylabel('Imaginary part [m/N]')
    plot.xlabel('Frequency [Hz]')
    plot.plot(f,imaginary_part,colors[i],label='\u03B6 = '+str(damping_ratio[i]))  # Imaginary part plot

    plot.subplot(2, 2, 4)
    plot.ylabel('Imaginary part')
    plot.xlabel('Real part')
    plot.plot(real_part,imaginary_part,colors[i],label='\u03B6 = '+str(damping_ratio[i]))  # Nyquist plot
plot.legend()
plot.show()
