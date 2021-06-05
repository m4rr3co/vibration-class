# Vibrations Course
# Mechanical Engineering Post Graduation Program
# UNESP - Ilha Solteira, SP - Brazil
# ####
# Work by Lucas Veronez Goulart Ferreira
# ####

import numpy as np
import matplotlib.pyplot as plot
from numpy.linalg import inv
from scipy.fft import fft,ifft

# System parameters
mass_structure = ms = 1  # Kilograms [kg]
mu = 0.1
mass_absorber = ma = mu*ms
stiffness_structure = ks = 1e4  # Newton/meter [N/m]
damping_ratio_structure = zetas = 0.001
damping_structure = cs = zetas*2*np.sqrt(ks*ms)
damping_ratio_opt = zetaa = np.sqrt((3/8)*(mu/((1+mu)**3)))
stiffness_absorber = stiffness_opt = ka = ((1/(1+mu))**2)*(ma*ks)/ms
damping_absorber = ca = zetaa*2*np.sqrt(ka*ma)

# Sampling parameters
sampling_interval = delta = 0.009
sampling_frequency = fs = 1/delta
tmax = 13
t = np.arange(0,tmax,delta)

# Frequency vector
f = np.linspace(0,fs,len(t))
ang_vel = w = 2*np.pi*f

# Calculating System Matrices
mass = m = np.array([[ms,0],[0,ma]])
damping = c = np.array([[cs+ca,-ca],[-ca,ca]])
stiffness = k = np.array([[ks+ka,-ka],[-ka,ka]])
Ha = []
Hs = []
for i in range(len(ang_vel)):
    impedance = Z = k-((w[i]**2)*m)+(1j*w[i])*c
    Z_inv = inv(Z)
    Hs.append(Z_inv[0][0])
    Ha.append(Z_inv[1][0])

plot.plot(f,np.abs(Ha))
plot.plot(f,np.abs(Hs))
plot.show()

double_sided_Ha = ds_Ha = Ha+np.conj(np.flip(Ha))  # Formation of the double-sided spectrum
double_sided_Hs = ds_Hs = Hs+np.conj(np.flip(Hs))  # Formation of the double-sided spectrum

# Input Random Force
Fs = np.random.normal(0,1,len(t))
Fs_w = fft(Fs)*delta

xa_t = ifft(ds_Ha*Fs_w)/delta
xs_t = ifft(ds_Hs*Fs_w)/delta

plot.plot(t,xa_t)
plot.plot(t,xs_t)
plot.show()
