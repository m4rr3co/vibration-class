# Vibrations Course
# Mechanical Engineering Post Graduation Program
# UNESP - Ilha Solteira, SP - Brazil
# ####
# Work by Lucas Veronez Goulart Ferreira
# ####
import numpy as np
from scipy.fft import fft
from scipy.signal import chirp

# Sampling parameters
sampling_interval = delta = 0.00072
sampling_frequency = fs = 1/delta
tmax = 331.1
t = np.arange(-tmax/2,tmax/2,delta)

# Signal parameters: Sine wave
f = 1  # Hz
x = np.sin(2*np.pi*f*t)
ft = fft(x)*delta

msv = np.mean(x**2)
int1 = sum(((np.abs(ft/tmax))**2))

# Difference between MSV of the signal and sum of of FFT's coefficients
print(str(msv-int1))  # -1.9324034702949966e-06

# Signal parameters: Chirp wave
x = chirp(t,f0=1,t1=tmax/2,f1=5)
ft = fft(x)*delta

msv = np.mean(x**2)
int1 = sum(((np.abs(ft)/tmax)**2))

# Difference between MSV of the signal and sum of of FFT's coefficients
print(str(msv-int1))  # -1.931049561976206e-06

# Signal parameters: Random signal
x = np.random.random(len(t))
ft = fft(x)*delta

msv = np.mean(x**2)
int1 = sum(((np.abs(ft)/tmax)**2))

# Difference between MSV of the signal and sum of of FFT's coefficients
print(str(msv-int1))  # -1.2883146801034862e-06
