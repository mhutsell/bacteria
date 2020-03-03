#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 14:39:43 2020

@authors: Mack, Camilla, Matthew
"""


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as sp
import scipy.integrate as odeint

data = pd.read_excel("../data/Bacteria_data.xlsx")

data.plot(kind = "scatter", x = "X", y="Y")
plt.figure()

y_vals = data["Y"]
x_vals = data["X"]
    
x_vals_sorted = sorted(x_vals)

y_vals_log = [np.log(y) for y in y_vals]




################### Section 1: Lin Fit for finding gmax and mu

x_vals_1 = x_vals[0:32]
y_vals_1 = y_vals_log[0:32]
np_y_1 = np.array(y_vals_1)
np_x_1 = np.array(x_vals_1)
np_x_1 = np_x_1.reshape((np_x_1.size, 1))
np_y_1 = np_y_1.reshape((np_y_1.size, 1))

a = float(np.linalg.lstsq(np_x_1, np_y_1, rcond=None)[0])
m, b = sp.linregress(x_vals_1, y_vals_1)[0:2]
print(m, b)

plt.plot(x_vals_1, y_vals_1, 'ro')
plt.figure()
y_vals_guessed_at = [a*x + b for x in x_vals_1]
plt.plot(x_vals_1, y_vals_guessed_at, 'ro')
print(a) # This equals gmax - mu, but mu is super small so probably just gmax
