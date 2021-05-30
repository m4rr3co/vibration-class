# Vibrations Course
# Mechanical Engineering Post Graduation Program
# UNESP - Ilha Solteira, SP - Brazil
# ####
# Work by Lucas Veronez Goulart Ferreira
# ####

import numpy as np
import matplotlib.pyplot as plot
from scipy import signal
from scipy.fft import fft

# System parameters
mass = m = 1  # Kilograms [kg]
stiffness = k = 1e4  # Newton/meter [N/m]
wn = np.sqrt(stiffness/mass)  # Natural frequency
critical_damping = cc = 2*m*wn  # zeta = c/cc
zeta = 0.001
c = zeta*cc
wd = wn*np.sqrt(1-(zeta ** 2))

# Sampling parameters
sampling_interval = delta = 0.079
sampling_frequency = fs = 1/delta
tmax = 133
t = np.arange(0,tmax,delta)

# Input Force
A = 100  # [N]
f_t = A*(np.random.rand(len(t))-0.5)
h = (np.exp(-zeta*wn*t)/wd*m)*np.sin(wd*t)
x_t = np.convolve(h,f_t,'full')
x_t_max = max(x_t)
x_t_noised = x_t[:len(t)] + 0.2*x_t_max*(np.random.rand(len(t))-0.5)

# Fourier Transforms
Xw = fft(x_t[:len(t)])*delta
Xw_noised = fft(x_t_noised)*delta
Fw = fft(f_t)*delta
Xw_double_sided = np.concatenate((np.flip(np.conj(Xw[:len(Xw)//2])),Xw[:len(Xw)//2+1]))
w = np.linspace(0,fs,len(t))
w_double_sided = np.linspace(-fs/2,fs/2,len(t))


# Power Spectral Densities
Sfx = (np.conj(Fw)*Xw)/tmax
Sff = (np.conj(Fw)*Fw)/tmax
Sxx = (np.conj(Xw)*Xw)/tmax

Sfx_noised = (np.conj(Fw)*Xw_noised)/tmax
Sff_noised = (np.conj(Fw)*Fw)/tmax
Sxx_noised = (np.conj(Xw_noised)*Xw_noised)/tmax

H1 = Sfx/Sff
H1_noised = Sfx_noised/Sff

gama_fx_squared = np.abs(Sfx**2)/(Sff*Sxx)
gama_fx_squared_noised = np.abs(Sfx_noised**2)/(Sff*Sxx_noised)

frf_analytical = 1 / (k-(m*(2*np.pi*w)**2)+((2*np.pi*w)*c*1j))

# plot.plot(w,np.abs(H1/tmax))
# plot.plot(w,np.abs(frf_analytical))
# plot.xlim([0,fs/2])
# plot.show()

# plot.plot(w,np.unwrap(np.angle(H1/tmax)))
# plot.plot(w,np.angle(frf_analytical))
# plot.xlim([0,fs/2])
# plot.show()

plot.plot(w,np.abs(gama_fx_squared))
plot.plot(w,np.abs(gama_fx_squared_noised))
plot.ylim([0.9,1.1])
plot.show()




