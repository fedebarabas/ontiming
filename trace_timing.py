# -*- coding: utf-8 -*-
"""
Created on Fri Mar 14 13:59:35 2014

@author: Federico Barabas

This is a simulation aiming to extract the time constant of a molecule's
intensity trace to get it's characteristic ontime.
"""

from __future__ import division, with_statement, print_function

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def autocorr(x):
    xc = np.correlate(x, x, mode='same') / (len(x) * np.mean(x)**2)
    return xc[xc.size/2:]


def exp_decay(x, A, lambd):
    return 1 + A * np.exp(- lambd * x)

# For now, fixed number of transitions and molecules
n_transitions = 10
n_molecules = 10

# Characteristic ontime and offtime of the molecule
ontime = 10      # [ms]
offtime = 1000   # [ms]
t_frame = 5      # [ms]

n = 200
koff = np.zeros(n)
kon = np.zeros(n)
x = np.arange(0, 25 * t_frame, t_frame)

for j in np.arange(n):

    # A single trace made of the concatenation of the traces of all molecules
    ontimes = np.random.exponential(ontime, n_transitions * n_molecules)
    offtimes = np.random.exponential(offtime, n_transitions * n_molecules)

    # Total time of the experiment
    n_frames = int(np.ceil(np.sum(ontimes + offtimes) / 100) * 100)

    # Trace calculation
    trace_real = np.zeros(100 * n_frames)
    t_transitions = np.cumsum(np.array(list(zip(ontimes, offtimes))).reshape(-1))

    trace_real[0:np.int(100 * t_transitions[0])] = 1
    for i in np.arange(1, len(t_transitions) / 2).astype(int):
        trace_real[np.int(100 * t_transitions[2 * i - 1]):
                   np.int(100 * t_transitions[2 * i])] = 1

    # Trace integration in each frame
    trace = np.array([np.sum(trace_real[100 * i:100 * (i + 1)])
                     for i in np.arange(0, n_frames, t_frame)]) / 100

    # Autocorrelation
    corr = autocorr(trace)
    c0 = corr[0]

    # Autocorrelation fitting
    k = (0.69 - np.log(1 - 1 / (c0 - 1))) / np.where(corr < c0 / 2)[0][0]
    guess = [corr[corr.argmax()] - 1, k]
    fit_par, fit_var = curve_fit(exp_decay, x, corr[:25], p0=guess)

    # kon and koff saving
    koff[j] = fit_par[1] * fit_par[0] / (fit_par[0] + 1)
    kon[j] = fit_par[1] * (1 - fit_par[0] / (fit_par[0] + 1))

time = np.arange(0, n_frames, t_frame)

# Plots
f, axarr = plt.subplots(2, 2)

# Measured trace
axarr[0, 0].plot(time, trace)
axarr[0, 0].grid('on')
axarr[0, 0].set_xlabel('Time [frames]')
axarr[0, 0].set_ylim([-0.2, 1.2])
axarr[0, 0].set_ylabel('On fraction')

# Autocorrelation of measured trace
axarr[1, 0].plot(time[:25], corr[:25], 'bo')
axarr[1, 0].plot(time[:25], exp_decay(x, fit_par[0], fit_par[1]), 'r')
axarr[1, 0].set_xlabel('Time [frames]')
axarr[1, 0].grid('on')

# Histogram of koff
axarr[0, 1].set_title('koff')
axarr[0, 1].hist(koff, bins=20)
axarr[0, 1].set_xlabel('Rate [frames^-1]')
axarr[0, 1].grid('on')

# Histogram of kon
axarr[1, 1].set_title('kon')
axarr[1, 1].hist(kon, bins=20)
axarr[1, 1].set_xlabel('Rate [frames^-1]')
axarr[1, 1].grid('on')
