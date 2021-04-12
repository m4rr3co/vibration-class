import numpy as np
from scipy import signal
import matplotlib.pyplot as plot
from scipy.fft import fft, fftfreq

amplitude = A = 1
period = 1
teta = np.arange(-3*period,3 * period,0.01)
square_wave = A * signal.square(2 * np.pi * (teta / period))

t1 = np.linspace(-3*period,3*period,len(square_wave),endpoint=False)
# plot.plot(t1,square_wave)
fourier = fft(square_wave)
freq = fftfreq(len(square_wave),0.01)
plot.plot(freq,np.abs(fourier))
plot.show()
