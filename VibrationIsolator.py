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

# System parameters
mass = m = 10  # Kilograms [kg]
fn = 5  # Natural Frequency [Hz]
wn = 2*np.pi*fn  # Natural Frequency [rad/s]
stiffness = k = (wn**2)*m  # Newton/meter [N/m]
print(k)
critical_damping = cc = 2*m*wn  # zeta = c/cc
zeta = 0.005
c = zeta*cc
wd = wn*np.sqrt(1-(zeta ** 2))

# Sampling parameters
sampling_interval = delta = 0.009
sampling_frequency = fs = 1/delta
tmax = 13
t = np.arange(0,tmax,delta)

# Frequency vector
f = np.linspace(0,fs,len(t))

# Displacement Input
x_max = 0.001  # [m]
x_t = x_max * (np.random.rand(len(t))-0.5)
Xw = fft(x_t)*delta

# Analytical Transfer Function
# Base Displacement Case
Hw = (k+1j*(2*np.pi*f)*c)/(k-((2*np.pi*f)**2)*m+(1j*(2*np.pi*f)*c))
Yw = Hw*Xw
y_t = ifft(Yw)/delta

# H estimate
H_est = Yw/Xw

# Extracting k
H_max = max(np.abs(H_est))
fn_est = f[np.where(np.abs(H_est) == H_max)]
print(m*(fn_est[0]*2*np.pi)**2)

# Extracting zeta
zeta_est = 1/(2*(H_max-1))
print(zeta)
print(zeta_est)
