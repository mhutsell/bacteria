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

################### Section 1: Lin Fit for finding mu
# We know that at small values of row, dn/dt = n(gmax-mu)

x_vals_1 = x_vals[40:]
y_vals_1 = y_vals_log[40:]

# np_y_1 = np.array(y_vals_1)
# np_x_1 = np.array(x_vals_1)
# np_x_1 = np_x_1.reshape((np_x_1.size, 1))
# np_y_1 = np_y_1.reshape((np_y_1.size, 1))

#Either works
#a = float(np.linalg.lstsq(np_x_1, np_y_1, rcond=None)[0])
m1, b1 = sp.linregress(x_vals_1, y_vals_1)[0:2]

print(m1, b1) # m and b in mx + b for our linear regression

plt.plot(x_vals_1, y_vals_1, 'ro')
plt.figure()
y_vals_guessed_at_1 = [m1*x + b1 for x in x_vals_1]
plt.plot(x_vals_1, y_vals_guessed_at_1, 'ro')
m1 *= -1

################### Section 2: Lin Fit for finding gmax - mu
# We know that at large values of row, dn/dt = n(gmax-mu)
x_vals_2 = x_vals[0:32]
y_vals_2 = y_vals_log[0:32]
np_y_1 = np.array(y_vals_1)

# np_x_1 = np.array(x_vals_1)
# np_x_1 = np_x_1.reshape((np_x_1.size, 1))
# np_y_1 = np_y_1.reshape((np_y_1.size, 1))

#Either works
#a = float(np.linalg.lstsq(np_x_1, np_y_1, rcond=None)[0])
m2, b2 = sp.linregress(x_vals_2, y_vals_2)[0:2]

print(m2, b2) # m and b in mx + b for our linear regression

plt.plot(x_vals_2, y_vals_2, 'ro')
plt.figure()
y_vals_guessed_at_2 = [m2*x + b2 for x in x_vals_2]
plt.plot(x_vals_2, y_vals_guessed_at_2, 'ro')

g_max_guess = m2+m1
print(g_max_guess) 

