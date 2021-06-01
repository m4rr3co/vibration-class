# Vibrations Course
# Mechanical Engineering Post Graduation Program
# UNESP - Ilha Solteira, SP - Brazil
# ####
# Work by Lucas Veronez Goulart Ferreira
# ####

import numpy as np
import matplotlib.pyplot as plot
from scipy.signal import hilbert

# Estimation of Damping from IRF
# IRF: x(t) = (e^(zeta*wn*t)/m*wd)*sin(wd*t)
# Parameters
mass = m = 1  # Kilograms [kg]
stiffness = 1e4  # Newton/meter [N/m]
wn = np.sqrt(stiffness/mass)  # Natural frequency

# Damping Ratio
# For now, these are given, but can be calculated from
# response as will be shown further on.
damping_ratio = [0.1,0.01,0.001]  # No dimension, i.e.: zeta = c/sqrt(m.k)
tmax = 2
delta = 0.001
interval = t = np.arange(0,tmax,delta)  # Time interval [0,0.001,...,10] seconds
upper_limit = 0.15
lower_limit = 0
secondary_interval = t2 = np.arange(lower_limit,upper_limit,0.001)

# Calculated Parameters
wd = [0,0,0]  # Defining array of values for Damped natural frequencies
for i in range(len(damping_ratio)):
    wd[i] = (wn*np.sqrt(1-(damping_ratio[i]**2)))
print('Damped natural frequencies [rad/s]: '+str(wd))

calculated_damping_ratios = cdr = [0,0,0]  # Defining array of values for Calculated Damping Ratios
# Calculated Damping Ratio: Angular Coefficient from Fitting Curve = -zeta*wn
# CDR = zeta = - ang.coef. / wn

colors = c = ['r','g','b']
lower_limits = [0.05,0.55,0.6]
upper_limits = [0.36,1.25,1.22]
for i in range(0,3):
    h = (np.exp(-damping_ratio[i]*wn*t)/wd[i]*m)*np.sin(wd[i]*t)
    signal_envelope = np.abs(hilbert(h))  # Absolute of Hilbert Transform of the signal

    # Chosen range for calculating the estimate damping value
    fit_range = signal_envelope[np.logical_and(t > lower_limits[i], t <= upper_limits[i])]
    range = np.arange(lower_limits[i],upper_limits[i],0.001)

    # Using polyfit in the specified range
    poly_fit = np.polyfit(range, np.log(fit_range), 1)
    expo_curve = np.exp(poly_fit[1]) * np.exp(poly_fit[0] * range)
    cdr[i] = -poly_fit[0]/wn
    plot.plot(t,signal_envelope,c[i],label='\u03B6 = '+str(damping_ratio[i]))
    plot.plot(np.arange(lower_limits[i],upper_limits[i],0.001),expo_curve)
plot.yscale('log')
plot.ylabel('Displacement [m]')
plot.xlabel('Time [s]')
plot.legend()
plot.show()
print('Real damping ratios (Zeta): '+str(damping_ratio))
print('Calculated damping ratios: '+str(cdr))
