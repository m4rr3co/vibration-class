# Vibrations Course
# Mechanical Engineering Post Graduation Program
# UNESP - Ilha Solteira, SP - Brazil
# ####
# Work by Lucas Veronez Goulart Ferreira
# ####

import numpy as np
import matplotlib.pyplot as plot
from scipy import signal

# Time domain convolution: f(t)*h(t) = x(t)
# 1st case: f(t) = delayed impulse

# System Parameters
mass = m = 1  # [Kg]
stiffness = k = 1e4  # [N/m]
wn = np.sqrt(k/m)  # [rad/s] - Natural Frequency
fn = wn / (2 * np.pi)  # [Hz] - Natural Frequency
damping_ratio = 0.01  # Zeta
critical_damping = cc = 2*np.sqrt(k*m)  # [Ns/m]
c = damping_ratio * critical_damping  # [Ns/m] - Damping
wd = wn * np.sqrt(1-damping_ratio**2)  # [rad/s] - Damped Natural Frequency
fd = wd / (2 * np.pi)  # [Hz] - Damped Natural Frequency
print('Calculated Damped Natural Frequency = '+str(fd)+' Hz.')

# Sampling Parameters
delta = 0.001
t_max = 10
time_range = t = np.arange(0,t_max,delta)

# 1st Case: Delayed impulse
h = (np.exp(-damping_ratio*wn*t)/wd*m)*np.sin(wd*t)  # Transfer Function
imp = signal.unit_impulse(len(time_range),int(len(time_range)*0.1))  # Delayed impulse

conv = np.convolve(h,imp,'full')
plot.subplot(4,2,1)
plot.plot(t,imp)
plot.ylabel('Unit Impulse')
plot.subplot(4,2,2)
plot.plot(t,conv[:len(t)])

# 2nd Case: Sine impulse
sine_imp = np.zeros(len(t))  # Padding vector with zeros
Tn = 2*np.pi/wn  # Natural Period
t_imp = np.arange(0,Tn,delta)
f_sine = 1/(2*Tn)
sine_signal = np.sin(2*np.pi*f_sine*t_imp)
for i in range(len(t_imp)):
    sine_imp[i] = sine_signal[i]

conv = np.convolve(h,sine_imp,'full')
plot.subplot(4,2,3)
plot.plot(t,sine_imp)
plot.ylabel('Sine Impulse')
plot.subplot(4,2,4)
plot.plot(t,conv[:len(t)])

# plot.plot(t,conv[:len(t)])
# plot.show()

# 3rd Case: Random input
rand_signal = 10*(np.random.rand(len(t))-0.5)  # Random input ranging [-5,5] N

conv = np.convolve(h,rand_signal,'full')
plot.subplot(4,2,5)
plot.plot(t,rand_signal)
plot.ylabel('Random Input')
plot.subplot(4,2,6)
plot.plot(t,conv[:len(t)])

# 4th Case: Chirp Signal
chirp_signal = signal.chirp(t,0,100,max(t))  # Random input with 5 N of amplitude

conv = np.convolve(h,chirp_signal,'full')
plot.subplot(4,2,7)
plot.plot(t,chirp_signal)
plot.ylabel('Chirp Signal')
plot.subplot(4,2,8)
plot.plot(t,conv[:len(t)])
plot.show()
