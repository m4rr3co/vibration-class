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

plot.plot(f,np.abs(Ha),label='Absorber')
plot.plot(f,np.abs(Hs),label='Structure')
plot.xlabel('Frequency [Hz]')
plot.ylabel('Xi/F \n i = Absorber, Structure')
plot.legend()
plot.show()

double_sided_Ha = ds_Ha = Ha+np.conj(np.flip(Ha))  # Formation of the double-sided spectrum
double_sided_Hs = ds_Hs = Hs+np.conj(np.flip(Hs))  # Formation of the double-sided spectrum

# Input Random Force
Fs = np.random.normal(0,1,len(t))
Fs_w = fft(Fs)*delta

xa_t = ifft(ds_Ha*Fs_w)/delta
xs_t = ifft(ds_Hs*Fs_w)/delta

# Fourier Transforms
Xaw = fft(xa_t)*delta
Xsw = fft(xs_t)*delta
Fw = Fs_w

# Power Spectral Densities Absorber
Sfxa = (np.conj(Fw)*Xaw)/tmax
Sffa = (np.conj(Fw)*Fw)/tmax
Sxxa = (np.conj(Xaw)*Xaw)/tmax

# Power Spectral Densities Structure
Sfxs = (np.conj(Fw)*Xsw)/tmax
Sffs = (np.conj(Fw)*Fw)/tmax
Sxxs = (np.conj(Xsw)*Xsw)/tmax

Ha1 = (Sfxa/Sffa)
Hs1 = (Sfxs/Sffs)

plot.subplot(2,2,1)
plot.plot(f,np.abs(ds_Ha),label='Absorber')
plot.plot(f,np.abs(Ha1),label='Structure')
plot.legend()
plot.ylabel('X/F')

plot.subplot(2,2,2)
plot.plot(f,np.abs(ds_Hs))
plot.plot(f,np.abs(Hs1))

plot.subplot(2,2,3)
plot.plot(f,np.unwrap(np.angle(ds_Ha)))
plot.plot(f,np.unwrap(np.angle(Ha1)))
plot.ylabel('Phase [rad]')
plot.xlabel('Frequency [Hz]')

plot.subplot(2,2,4)
plot.plot(f,np.unwrap(np.angle(ds_Hs)))
plot.plot(f,np.unwrap(np.angle(Hs1)))
plot.xlabel('Frequency [Hz]')

plot.show()
