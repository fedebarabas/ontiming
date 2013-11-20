# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 14:49:49 2013

@author: fbaraba
"""

from __future__ import division, with_statement, print_function

import numpy as np
#import matplotlib.pyplot as plt

from pint import UnitRegistry
ureg = UnitRegistry()
Q_ = ureg.Quantity

def dgauss(x, y, center, fwhm):
    return np.exp(-(x - center)**2 / (2 * fwhm**2)) * np.exp(-(y - center)**2 / (2 * fwhm**2))

# Constant definitions
nm_per_pixel = 133.0
ureg.define('pixel = {} * nanometer = px'.format(nm_per_pixel))

# Parameters
trc_window_l = 5 * ureg.px
psf_width = 250 * ureg.nanometer
emission_rate = 1000 / ureg.millisecond     # Means 1000 photons per ms
ontime = 5 * ureg.millisecond
exposure = 1 * ureg.millisecond

# Pixel grid
x = np.arange(0, trc_window_l.magnitude)
y = x
xx, yy = np.meshgrid(x, y, sparse=True)

# PSF calculation
psf = 3 * dgauss(xx, yy,
                   np.floor(trc_window_l.magnitude / 2),
                   psf_width.to(ureg.px).magnitude)

# Time goes from 0 ms to 3 * ontime with a step two orders of
# magnitude better than the exposure time
resolution = exposure / 100
t_window = 3 * ontime
time = np.arange(0, t_window.magnitude, resolution.magnitude)

t_ons = np.arange(4, 5, 0.1) * ureg.millisecond

n = 100
ratio = np.zeros((t_ons.size, n))

# Windows with gaussian noise incorporated
windows = np.random.normal(loc=0,
                           scale=5,
                           size=(t_ons.size,
                                 n,
                                 time.size,
                                 trc_window_l.magnitude,
                                 trc_window_l.magnitude))

# Photon's contribution to signal
win_photons = np.random.poisson(lam=psf,
                                size=(t_ons.size,
                                      n,
                                      (ontime / resolution).magnitude,
                                      trc_window_l.magnitude,
                                      trc_window_l.magnitude))

for i in range(len(t_ons)):
    init = round((t_ons[i] / resolution).magnitude)
    end = round(((t_ons[i] + ontime)/ resolution).magnitude)
    # noise + signal when dye is on
    windows[i, :, init:end] = windows[i, :, init:end] + win_photons[i]
    
# Reshape grouping signal from each exposure
windows_ccd = windows.reshape(t_ons.size, 
                              n, 
                              t_window.magnitude, 
                              -1, 
                              x.size, x.size)
                              
signal = np.round(windows_ccd.sum(3))

# CCD plot
#h = plt.imshow(signal[7], cmap=plt.cm.jet, interpolation='nearest')
#cbar = plt.colorbar()
#plt.show()

trace = signal.sum(axis=(3,4))
ratio = trace[:, :, 4] / trace[:, :, 5]



for i in range(len(ratio)):
    print(ratio[i].mean())
    print(ratio[i].std())
#    plt.plot(trace[i, 0])

#plt.grid('on')
#plt.xlabel('Time [{}]'.format(ontime.units))
#plt.ylabel('Photons')
#plt.show()