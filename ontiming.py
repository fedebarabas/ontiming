# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 14:49:49 2013

@author: fbaraba
"""

from __future__ import division, with_statement, print_function

import numpy as np
import matplotlib.pyplot as plt

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
psf = 250 * dgauss(xx, yy,
                   np.floor(trc_window_l.magnitude / 2),
                   psf_width.to(ureg.px).magnitude)

# Time goes from 0 ms to 3 * ontime with a step two orders of
# magnitude better than the exposure time
resolution = exposure / 100
time = np.arange(0, 3 * ontime.magnitude, resolution.magnitude)
state = np.zeros(time.shape)
t_on = 5 * ureg.millisecond
state[(t_on / resolution).magnitude:
     ((t_on + ontime)/ resolution).magnitude] = 1

# Adding noise
readout = np.random.normal(loc=50, scale=10, size=psf.shape)
window = np.array([np.random.poisson(lam=psfi) for psfi in psf])
#window = readout
window = window + readout

# Plot
h = plt.imshow(window, cmap=plt.cm.jet, interpolation='nearest')
cbar = plt.colorbar()
plt.show()

print(np.sum(window))