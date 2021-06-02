# Vibrations Course
# Mechanical Engineering Post Graduation Program
# UNESP - Ilha Solteira, SP - Brazil
# ####
# Work by Lucas Veronez Goulart Ferreira
# ####

import numpy as np
import matplotlib.pyplot as plot
from scipy import signal
from scipy.fft import fft,ifft

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
sampling_frequency = fs = 1/delta
t_max = 10
time_range = t = np.arange(0,t_max,delta)

freq_bins = f = np.linspace(0,fs,len(t))

# 1st Case: Delayed impulse
h = (np.exp(-damping_ratio*wn*t)/wd*m)*np.sin(wd*t)  # Transfer Function
imp = signal.unit_impulse(len(time_range),int(len(time_range)*0.1))  # Delayed impulse

fft_h = fft(h)
fft_imp = fft(imp)

output = ifft(fft_h*fft_imp)

plot.subplot(4,1,1)
plot.ylabel('Unit Impulse')
plot.plot(t,np.real(output))

# 2nd Case: Sine impulse
sine_imp = np.zeros(len(t))  # Padding vector with zeros
Tn = 2*np.pi/wn  # Natural Period
t_imp = np.arange(0,Tn,delta)
f_sine = 1/(2*Tn)
sine_signal = np.sin(2*np.pi*f_sine*t_imp)
for i in range(len(t_imp)):
    sine_imp[i] = sine_signal[i]

fft_h = fft(h)
fft_imp = fft(sine_imp)

output = ifft(fft_h*fft_imp)

plot.subplot(4,1,2)
plot.ylabel('Sine Impulse')
plot.plot(t,np.real(output))

# plot.plot(t,conv[:len(t)])
# plot.show()

# 3rd Case: Random input
rand_signal = 10*(np.random.rand(len(t))-0.5)  # Random input ranging [-5,5] N

fft_h = fft(h)
fft_imp = fft(rand_signal)

output = ifft(fft_h*fft_imp)

plot.subplot(4,1,3)
plot.ylabel('Random Input')
plot.plot(t,np.real(output))

# 4th Case: Chirp Signal
chirp_signal = signal.chirp(t,0,100,max(t))  # Random input with 5 N of amplitude

fft_h = fft(h)
fft_imp = fft(chirp_signal)

output = ifft(fft_h*fft_imp)

plot.subplot(4,1,4)
plot.ylabel('Chirp Signal')
plot.xlabel('Time [s]')
plot.plot(t,np.real(output))
plot.show()
