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

def uno(x, y):
    return 1

# Constant definitions
nm_per_pixel = 133.0
ureg.define('pixel = {} * nanometer = px'.format(nm_per_pixel))

trc_window_l = 5 * ureg.px
psf_width = 250 * ureg.nanometer
emission_rate = 1000 / ureg.millisecond     # Means 1000 photons per ms

# PSF calculation
x = np.arange(0, trc_window_l.magnitude, 0.1)
y = np.arange(0, trc_window_l.magnitude, 0.1)
xx, yy = np.meshgrid(x, y, sparse=True)
psf = dgauss(xx, yy, trc_window_l.magnitude / 2, psf_width.to(ureg.px).magnitude)

# Plot
A = psf.reshape(5,5,-1)
A.sum(-1)
print(A)
#h = plt.contourf(x, y, psf)
h = plt.contourf(np.arange(0, 5, 1), np.arange(0, 5, 1), A.sum(-1))
##h = plt.imshow(psf[25:25, 25:25], cmap=plt.cm.gray, interpolation='nearest')
h = plt.imshow(A.sum(-1), cmap=plt.cm.gray, interpolation='nearest')
cbar = plt.colorbar()
plt.grid()




print(psf)

# This is no good, I have to change to theoretical PSF
#counts_x = np.random.normal(loc=trc_window_l.magnitude / 2,
#                            scale=psf_width.to(ureg.px).magnitude,
#                            size=(5000))
#counts_y = np.random.normal(loc=trc_window_l.magnitude / 2,
#                            scale=psf_width.to(ureg.px).magnitude,
#                            size=(5000))
#
##plt.scatter(counts_x, counts_y)
#H, xedges, yedges = np.histogram2d(counts_x, counts_y,
#                                   range=[[0, 5], [0, 5]], bins=5)
#extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]
#plt.imshow(H, extent=extent, interpolation='None')
#plt.colorbar()