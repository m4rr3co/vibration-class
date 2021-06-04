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
from scipy.signal import coherence as ch

# System parameters
mass = m = 1  # Kilograms [kg]
stiffness = k = 1e4  # Newton/meter [N/m]
wn = np.sqrt(stiffness/mass)  # Natural frequency
critical_damping = cc = 2*m*wn  # zeta = c/cc
zeta = [0.1,0.01,0.001]
c = [0,0,0]
for i in range(len(zeta)):
    c[i] = zeta[i]*cc
wd = [0,0,0]
for i in range(len(zeta)):
    wd[i] = wn*np.sqrt(1-(zeta[i] ** 2))

# Sampling parameters
sampling_interval = delta = 0.001
sampling_frequency = fs = 1/delta
tmax = 5
t = np.arange(0,tmax,delta)

for i in range(len(zeta)):
    # Input Force
    f_t = np.random.normal(0,1,len(t))
    h = (np.exp(-zeta[i]*wn*t)/wd[i]*m)*np.sin(wd[i]*t)
    x_t = np.convolve(h,f_t,'full')[:len(t)]
    x_t = x_t + np.random.normal(0,0.00005,len(t))
    noise = max(x_t)

    # Fourier Transforms
    freq_bins = f = np.linspace(0,fs,len(t))
    Xw = fft(x_t)*delta
    Fw = fft(f_t)*delta

    # Power Spectral Densities
    Sfx = (np.conj(Fw)*Xw)/tmax
    Sff = (np.conj(Fw)*Fw)/tmax
    Sxx = (np.conj(Xw)*Xw)/tmax

    H1 = Sfx/Sff

    f1,cohe = ch(f_t,x_t,1/delta,nperseg=2048)
    plot.plot(f1,cohe)
    plot.ylim([0,1.1])

    # frf_analytical = 1 / (k-(m*(2*np.pi*f)**2)+((2*np.pi*f)*c*1j))

    # plot.plot(w,np.abs(H1/tmax))
    # plot.plot(w,np.abs(frf_analytical))
    # plot.xlim([0,fs/2])
    # plot.show()

    # plot.plot(w,np.unwrap(np.angle(H1/tmax)))
    # plot.plot(w,np.angle(frf_analytical))
    # plot.xlim([0,fs/2])
    # plot.show()
plot.show()




