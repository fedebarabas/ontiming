# -*- coding: utf-8 -*-
"""
Created on Fri Mar 14 13:59:35 2014

@author: Federico Barabas

This is a simulation aiming to extract the time constant of a molecule's
intensity trace to get it's characteristic ontime.
"""

from __future__ import division, with_statement, print_function

import numpy as np
import matplotlib.pyplot as plt


def autocorr(x):
    result = np.correlate(x, x, mode='same')
    return result[result.size/2:]

# For now, fixed number of transitions and molecules
n_transitions = 10
n_molecules = 15

# Characteristic ontime and offtime of the molecule, in frame units
ontime = 5
offtime = 500

# A single trace made of the concatenation of the traces of all molecules
ontimes = np.random.exponential(ontime, n_transitions * n_molecules)
#ontimes = ontime * np.ones(n_transitions * n_molecules)
#ontimes = np.random.normal(ontime, ontime, n_transitions * n_molecules)
#offtimes = np.random.normal(offtime, offtime, n_transitions * n_molecules)
#offtimes = offtime * np.ones(n_transitions * n_molecules)
offtimes = np.random.exponential(offtime, n_transitions * n_molecules)

# Total time of the experiment
t_total = np.ceil((ontimes.sum() + offtimes.sum()) / 100) * 100

# Trace calculation
trace_real = np.zeros(100 * t_total)
trace = np.zeros(t_total)
t_transitions = np.cumsum(np.array(zip(ontimes, offtimes)).reshape(-1))
trace_real[0:100 * t_transitions[0]] = 1
for i in np.arange(1, len(t_transitions) / 2):
    trace_real[100 * t_transitions[2 * i - 1]:100 * t_transitions[2 * i]] = 1

# Trace integration in each frame
trace = np.array([np.sum(trace_real[100 * i:100 * (i + 1)])
                 for i in np.arange(t_total)]) / 100

# Autocorrelation
corr = autocorr(trace)

# Plots
f, axarr = plt.subplots(2, 1)

# Emission trace
#axarr[0].plot(np.arange(0, t_total, 0.01), trace_real)
#axarr[0].grid('on')
#axarr[0].set_xlabel('Time [frames]')
#axarr[0].set_ylim([-0.2, 1.2])
#axarr[0].set_ylabel('State')

# Measured trace
axarr[0].plot(trace)
axarr[0].grid('on')
axarr[0].set_xlabel('Time [frames]')
axarr[0].set_ylim([-0.2, 1.2])
axarr[0].set_ylabel('On fraction')

# Autocorrelation of measured trace
axarr[1].plot(corr)
axarr[1].set_xlabel('Time [frames]')
axarr[1].grid('on')
