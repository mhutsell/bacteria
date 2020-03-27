#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 14:39:43 2020

@authors: Mack
"""

import pprint
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as sp
from scipy.integrate import odeint
from scipy.optimize import curve_fit

data = pd.read_excel("../data/Bacteria_data.xlsx")
data.plot(kind = "scatter", x = "X", y="Y")
plt.figure()

y_vals = data["Y"]
x_vals = data["X"]
    
x_vals_sorted = sorted(x_vals)

y_vals_log = [np.log(y) for y in y_vals]

################### Section 1: Lin Fit for finding mu
# We know that at small values of row, dn/dt = n(mu)

x_vals_1 = x_vals[40:]
y_vals_1 = y_vals_log[40:]

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

m2, b2 = sp.linregress(x_vals_2, y_vals_2)[0:2]

print(m2, b2) # m and b in mx + b for our linear regression

plt.plot(x_vals_2, y_vals_2, 'ro')
plt.figure()
y_vals_guessed_at_2 = [m2*x + b2 for x in x_vals_2]
plt.plot(x_vals_2, y_vals_guessed_at_2, 'ro')

g_max_guess = m2+m1
print(g_max_guess) 

#################### Section 3: Finding a
# We know that, disregarding death of bacteria as we're just approximating, our
# amount of nutrient, 0.2 microgram/ml, led to the creation of roughly 2e8
# bacteria. Therefore, it takes 0,2 microgram/ml to make 2e8 bacteria per ml.
# so we can approximate 'a' by taking 2e8 and dividing it by 0.2. This will
# tell us how many bacteria one microgram/ml creates.

a_guess = 0.2 / np.max(y_vals)
print(a_guess)

#################### Section 4: Approximating nutrient decay


INIT_NUTR_AMOUNT = 0.2          # microgram / ml
nutr_amount = INIT_NUTR_AMOUNT


######## Finding the average y value for each x value
x_val_avg_y = {}
for x in x_vals:
    x_val_avg_y[x] = 0
for i, x in enumerate(x_vals):
    x_val_avg_y[x] += y_vals[i]

for x in set(x_vals):
    x_val_avg_y[x] /= len([y for y in x_vals if y == x])

# =============================================================================
# #### Finding the average y value for each x value
# init_ind = 0
# x_val_avg_y = {}
# for i in range(1,52):
#     total_bact = 0
#     moving_ind = i-1
#     if not (x_vals[moving_ind] == x_vals[i]):
#         while (moving_ind >= init_ind):
#             total_bact += y_vals[moving_ind]
#             moving_ind -= 1
#         x_val_avg_y[x_vals[i-1]] = float(total_bact) / (i-1 - init_ind)
# 
# x_val_avg_y[100] = (float(y_vals[50]) + y_vals[51]) / 2
# =============================================================================

# Now approximating nutrient decay
nutrient_concentration_x = set(x_vals)
nutrient_dict = {}
for x in nutrient_concentration_x:
    if x > 40:
        nutrient_dict[x] = 0.001
    else:
        nutrient_dict[x] = 0.2 - float(x_val_avg_y[x]) * a_guess


approx_row_d = [(y_vals[i+1]- y_vals[i]) for i in range(0,51)]
approx_row_d +=[approx_row_d[50]] # so that we can have 52 values
k_guess_list = []
for i in range(32,41):   
    k_guess_at_i = ((y_vals[i] * g_max_guess * nutrient_dict[x_vals[i]]) / (y_vals[i] * m1 + approx_row_d[i]) - nutrient_dict[x_vals[i]]) 
    k_guess_list.append(k_guess_at_i)

k_guess = np.average(k_guess_list)


# other nutrient dict w/ full values:
nut_list_2 = []
for i,x in enumerate(x_vals):
    nut_list_2.append(0.2 - y_vals[i] * a_guess)

#################### Section 5: curve_fit

init_y = [x_val_avg_y[x_vals[0]], 0.2]

def diff_n_p(x, t, gmax, k, mu, a):
    n, row = x
    dndp = [(n*gmax*row)/(row+k) - (n*mu), (-a)*(n*gmax*row)/(row+k)]  
    return dndp

def log_diff(t, g_max_guess, k_guess, m1, a_guess):
    result = odeint(diff_n_p, init_y, t, args=(g_max_guess, k_guess, m1, a_guess)).flatten()
    for i, x in enumerate(result):
        if x <= 0:
            result[i] = 0.001
    return np.log(result)


pop_log_np = np.log(np.array(list(x_val_avg_y.values())))
row_log_np = np.log(np.array(list(nutrient_dict.values())))
y_vals = np.array([pop_log_np, row_log_np]).transpose()

x_np = np.array(list(x_val_avg_y.keys()))

g_max, k, mu, a = curve_fit(log_diff, x_np, y_vals.flatten(), p0=(g_max_guess, k_guess, m1, a_guess))[0]
print(g_max, k, mu, a)
print(g_max_guess, k_guess, m1, a_guess)


t = np.arange(4.5,150.1,0.1)
def modeling(x, t):
    n, row = x
    dndp = [(n*g_max*row)/(row+k) - (n*mu), (-a)*(n*g_max*row)/(row+k)]
    return dndp

result_1 = odeint(modeling, init_y, t)
