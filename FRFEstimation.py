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
    x_t_noised = x_t + np.random.normal(0,0.00005,len(t))

    # Fourier Transforms
    freq_bins = f = np.linspace(0,fs,len(t))
    Xw = fft(x_t)*delta
    Fw = fft(f_t)*delta

    # Power Spectral Densities
    Sfx = (np.conj(Fw)*Xw)/tmax
    Sff = (np.conj(Fw)*Fw)/tmax
    Sxx = (np.conj(Xw)*Xw)/tmax

    # H1 = (Sfx/Sff)*delta
    # frf_analytical = fft(h)*delta
    # plot.subplot(2,1,1)
    # plot.title('Transfer Function for \u03B6 = 0.1')
    # plot.plot(f, np.abs(H1),label='TF Estimate')
    # plot.plot(f,np.abs(frf_analytical),label='FFT of h(t)')
    # plot.ylabel('Transfer Function H = X/F')
    # plot.legend()
    # plot.subplot(2,1,2)
    # plot.plot(f,np.angle(H1))
    # plot.plot(f,np.unwrap(np.angle(frf_analytical)))
    # plot.xlabel('Frequency [Hz]')
    # plot.ylabel('Phase \u03C6 [rad]')
    # plot.show()

    f1,cohe = ch(f_t,x_t,1/delta,nperseg=2048)
    f2,cohe2 = ch(f_t, x_t_noised, 1 / delta, nperseg=2048)
    plot.plot(f2,cohe2,label='Noised')
    plot.plot(f1, cohe, label='Non-noised')
    plot.ylim([0,1.1])
    plot.ylabel('\u03B3xy^2')
    plot.xlabel('Frequency [Hz]')
    plot.legend()
    plot.show()
