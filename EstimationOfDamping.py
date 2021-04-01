# Vibrations Course
# Mechanical Engineering Post Graduation Program
# UNESP - Ilha Solteira, SP - Brazil
# ####
# Work by Lucas Veronez Goulart Ferreira
# ####

import numpy as np
import matplotlib.pyplot as plot
from scipy.signal import hilbert

# Estimation of Damping
# IRF: x(t) = (e^(zeta*wn*t)/m*wd)*sin(wd*t)
# Parameters
mass = m = 1  # Kilograms [kg]
stiffness = 10e4  # Newton/meter [N/m]
wn = np.sqrt(stiffness/mass)  # Natural frequency

# Damping Ratio
# For now, these are given, but can be calculated from
# response as will be shown further on.
damping_ratio = [0.1,0.01,0.001]  # No dimension, i.e.: zeta = c/sqrt(m.k)
pi = np.pi  # Pi = Everyone knows what this is
interval = t = np.arange(0.001,1.5,0.001)  # Time interval [0,0.001,...,10] seconds
upper_limit = 0.15
lower_limit = 0
secondary_interval = t2 = np.arange(lower_limit,upper_limit,0.001)

# Calculated Parameters
wd = [0,0,0]  # Defining array of values for Damped natural frequencies
for i in range(len(damping_ratio)):
    wd[i] = wn*np.sqrt(1-np.float_power(damping_ratio[i],2))
print('Damped natural frequencies [Hz]: '+str(wd))

calculated_damping_ratios = cdr = [0,0,0]  # Defining array of values for Calculated Damping Ratios
# Calculated Damping Ratio: Angular Coefficient from Fitting Curve = -zeta*wn
# CDR = zeta = - ang.coef. / wn

colors = c = ['r','g','b']
for i in range(0,3):
    x = (np.exp(-damping_ratio[i]*wn*t)/wd[i]*m)*np.sin(wd[i]*t)
    hilbert_result = hilbert(x)  # Hilbert Transform of the signal
    signal_envelope = np.abs(hilbert_result)  # Envelope of the signal
    fit_range = signal_envelope[np.logical_and(t > lower_limit, t <= upper_limit)]
    poly_fit = np.polyfit(t2, np.log(fit_range), 1)
    expo_curve = np.exp(poly_fit[1]) * np.exp(poly_fit[0] * t)
    cdr[i] = -poly_fit[0]/wn
    plot.plot(t,signal_envelope,c[i])
    plot.plot(t,expo_curve,c[i])
# plot.xscale('log')
plot.yscale('log')
plot.show()
print('Damping ratios (Zeta): '+str(damping_ratio))
print('Calculated damping ratios: '+str(cdr))





