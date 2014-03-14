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

# For now, fixed number of transitions and molecules
n_transitions = 10
n_molecules = 20

# Characteristic ontime and offtime of the molecule, in frame units
ontime = 3
offtime = 10

# A single trace made of the concatenation of the traces of all molecules
ontimes = np.random.exponential(ontime, n_transitions * n_molecules)
offtimes = np.random.exponential(offtime, n_transitions * n_molecules)

t_total = np.ceil((ontimes.sum() + offtimes.sum()) / 1000) * 1000

trace_real = np.array(zip(ontimes, offtimes)).reshape(-1)
trace = np.zeros(t_total)

ti = 0
for i in np.arange(t_total):

    t_on = 0
    t_off = 0

    while ti < i + 1:
        t_on +=
