# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 14:49:49 2013

@author: fbaraba
"""

from __future__ import division, with_statement, print_function

import numpy as np
import matplotlib.pyplot as plt
import airygauss as ag


def get_psf(x, y, em_rate, exposure, fwhm):
    center = x[0, 2]
    photons = em_rate * exposure
    A = 4 * np.log(2) * photons / (np.pi * fwhm**2)
    if type(A) is np.ndarray:
        return np.array([Ai * np.exp(-4 * np.log(2) *
                                     (x - center)**2 / fwhm**2) *
                        np.exp(-4 * np.log(2) *
                               (y - center)**2 / fwhm**2) for Ai in A])

    else:
        return A * np.exp(-4 * np.log(2) *
                          (x - center)**2 / fwhm**2) * np.exp(-4 * np.log(2) *
                                                              (y - center)**2
                                                              / fwhm**2)

# Constant definitions
nm_per_pixel = 133.0    # [nm]
em_wavelength = 670     # [nm]
NA = 1.4

# Parameters
trc_window_l = 5        # [px]
em_rate = 1000          # [1/ms]
ontime = 5              # [ms]
exposure = 1            # [ms]

#def ontiming():

# Pixel grid
x = np.arange(0, trc_window_l)
y = x
xx, yy = np.meshgrid(x, y, sparse=True)

# Time variables
t_window = 3 * ontime
time = np.arange(0, t_window, exposure)
t_ons = np.arange(4, 5, 0.05)        # [ms]
t_offs = t_ons + ontime

# p is the ratio of frame with dye emission in each case
p_on = np.transpose([time[np.where(time > t_on)[0][0]] - t_on
                    for t_on in t_ons])

p_off = np.transpose([t_off - time[np.where(time > t_off)[0][0] - 1]
                     for t_off in t_offs])

p_m = np.ones((ontime / exposure - 1, t_ons.size))
p = np.vstack((p_on, p_m, p_off)).transpose()

# PSF calculation
fwhm = ag.fwhm(em_wavelength, NA) / nm_per_pixel    # [px]
psfs = np.array([get_psf(xx, yy, em_rate, pi, fwhm) for pi in p])

trials = 10000

# Photon's contribution to signal
photons = np.random.poisson(lam=psfs, size=(trials,
                                            t_ons.size,
                                            p.shape[1],
                                            trc_window_l,
                                            trc_window_l))

# Windows with gaussian noise incorporated
boxes = np.round(np.random.normal(loc=0, scale=20, size=(trials,
                                                         t_ons.size,
                                                         time.size,
                                                         trc_window_l,
                                                         trc_window_l)))

# If you want it without gaussian noise
#windows = np.zeros((trials, t_ons.size, time.size, trc_window_l,
#                     trc_window_l))

boxes[:, :, t_ons[0]:t_ons[0] + ontime + 1] += photons

trace = boxes.sum((3, 4))

#
#    # CCD plot
#    #h = plt.imshow(signal[7], cmap=plt.cm.jet, interpolation='nearest')
#    #cbar = plt.colorbar()
#    #plt.show()
#
ratio = np.zeros((trials, t_ons.size))
ratio = trace[:, :, 4] / trace[:, :, 5:9].mean()

f, axarr = plt.subplots(2, 2)
axarr[0, 0].plot(t_ons, ratio.std(0))
axarr[0, 0].grid('on')
axarr[0, 0].set_xlabel('On time [ms]')
axarr[0, 0].set_ylabel('Precision in on time determination [ms]')


axarr[1, 0].errorbar(t_ons, 5 - ratio.mean(0), yerr=ratio.std(0), fmt='o')
axarr[1, 0].grid('on')
axarr[1, 0].set_xlabel('On time [ms]')
axarr[1, 0].set_ylabel('Predicted time [ms]')

for i in range(trace.shape[1]):
    axarr[0, 1].plot(trace[0, i])

axarr[0, 1].grid('on')
axarr[0, 1].set_xlabel('Time [ms]')
axarr[0, 1].set_ylabel('Photons')
plt.show()
