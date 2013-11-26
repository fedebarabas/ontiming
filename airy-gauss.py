# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 18:04:32 2013

@author: fbaraba
"""

from __future__ import division, with_statement, print_function

import numpy as np
from scipy.special import jn
from scipy.optimize import curve_fit

import matplotlib.pyplot as plt

def gauss(x, fwhm):
    return np.exp(- 4 * np.log(2) * (x / fwhm)**2)

def airy(x):
    return (2 * jn(1,2 * np.pi * x) / (2 * np.pi * x))**2

x = np.arange(-2, 2, 0.01)
y = airy(x)

# Fitting only inside first Airy's ring
fit_int = np.where(abs(x) < 0.61)[0]

fit_par, fit_var = curve_fit(gauss,
                             x[fit_int[0]:fit_int[-1]],
                             y[fit_int[0]:fit_int[-1]],
                             p0=0.5)
fit = gauss(x, fit_par)

print(fit_par[0], fit_var[0,0])

plt.plot(x, y)
plt.plot(x[fit_int[0]:fit_int[-1]], fit[fit_int[0]:fit_int[-1]])
plt.grid('on')
plt.show()