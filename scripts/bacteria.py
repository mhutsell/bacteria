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

data = pd.read_excel("../data/Bacteria_data.xlsx") # The file path here will have to be the same on your computer or
                                                   # the code will have an error. I formatted my excel sheet to have
                                                   # two columns— one with x at the top and one with y at the top.
                                                   # And I sorted the data in the excel sheet so that I wouldn't have
                                                   # to do it here. The excel sheet will be attached in our email.

data.plot(kind = "scatter", x = "X", y="Y") # Scatter plot our experimental data. This corresponds to Figure 2.

y_vals = data["Y"] # get a list of our 'y' values 
x_vals = data["X"] # and 'x' values separately. It makes the code a bit more readable.

y_vals_log = [np.log(y) for y in y_vals] # Take the log of all of our y values.

################### Section 1: Lin Fit for finding mu
# We know that at small values of ro, dn/dt approximately equals n(mu)
# Therefore, we can use linear regression to get an educated guess for mu.
# Since bacteria population is declining in later time values of our data, we can assume that ro is
# small. Not a perfect calculation, but good for getting a guess.

x_vals_1 = x_vals[40:] # We are taking data points 40-51, as they are when population is declining most clearly
y_vals_1 = y_vals_log[40:] # We take log Y values because we expect that population decline will be exponential.

m1, b1 = sp.linregress(x_vals_1, y_vals_1)[0:2] # Then we do linear regression to find m1 (our mu guess) and b1.

print(m1)

# Plot log of y values vs x values for data points 40-51.
plt.figure()
plt.plot(x_vals_1, y_vals_1, 'ro')
plt.title("Log Y Values vs X values 40-51")
plt.xlabel("Hours")
plt.ylabel("log(bacterial concentration) (CFU/ml)")

# Plot guessed y values vs x values for x values 40-51. Guessed y values are guesses of log(bacterial concentration).
# Guesses are made using m1 and b1 as previously calculated.
plt.figure()
y_vals_guessed_at_1 = [m1*x + b1 for x in x_vals_1] # Guess Y values 40-51 based on our mu_guess.
plt.plot(x_vals_1, y_vals_guessed_at_1, 'ro', color = "blue")
plt.title("Guessed Log Y Values vs X Values 40-51")
plt.xlabel("Hours")
plt.ylabel("log(bacterial concentration) (CFU/ml)")

m1 *= -1 # Our calculated m1 is negative, which is correct, but in our differential equation, we subtract
         # n(mu). Therefore, to avoid two negatives becoming a positive, I make m1 a positive value by multiplying
         # by negative 1. 

################### Section 2: Lin Fit for finding gmax - mu
# We know that at large values of ro, dn/dt = n(gmax-mu)
# At the beginning, ro is the largest that it will ever be, so the beginning has to have 'large' values of ro (for our data).
# To know when to stop taking data points, we simply looked at where the bacteria population began to stop displaying exponential
# growth. At that point, we can assume ro has become smaller and is no longer 'larger'. I put 'large' in quotes, because it being
# large is relative to our bacteria population's size, and so its exact numerical value can change from data set to data set.

# Take x and log(y) values from 0-31.
# We take log values of y because we assume growth to be exponential, and so this will allow us to fit linearly.  
x_vals_2 = x_vals[0:32]
y_vals_2 = y_vals_log[0:32]

m2, b2 = sp.linregress(x_vals_2, y_vals_2)[0:2] # Linear regression to find guesses for m2 and b2.
                                                # m2 is our guess for gmax-mu

print(m2) # m and b in mx + b for our linear regression

# 
plt.figure()
plt.plot(x_vals_2, y_vals_2, 'ro')
plt.title("Log Y Values vs X Values 0-31")
plt.xlabel("Hours")
plt.ylabel("log(bacterial concentration) (CFU/ml)")

plt.figure()
y_vals_guessed_at_2 = [m2*x + b2 for x in x_vals_2]
plt.plot(x_vals_2, y_vals_guessed_at_2, 'ro', color = "blue")
plt.title("Guessed Log Y Values vs X Values 0-31")
plt.xlabel("Hours")
plt.ylabel("log(bacterial concentration) (CFU/ml)")

g_max_guess = m2+m1 # m2 = gmax - mu, so to find our gmax guess, we must add our mu guess (m1)
print(g_max_guess) # Our gmax guess

#################### Section 3: Finding a
# We know that, disregarding death of bacteria as we're just approximating, our
# amount of nutrient, 0.2 microgram/ml, led to the creation of roughly 2e8
# bacteria. Therefore, it takes 0,2 microgram/ml to make 2e8 bacteria per ml.
# so we can approximate 'a' by taking 2e8 and dividing it by 0.2. This will
# tell us how many bacteria one microgram/ml creates.

a_guess = 0.2 / np.max(y_vals)
print(a_guess) # Our a guess

#################### Section 4: Approximating nutrient decay due to bacterial growth


######## Finding the average y value for each x value
# We do this because when construct our "remaining nutrient" dictionary,
# It made more sense to us to have an average y value for each x value,
# So that we didn't have wildly different amounts for remaining nutrients
# on the same x value. It also made using a dictionary possible, as individual
# keys were unique vs trying to have a dictionary with 2 '4' keys. 
x_val_avg_y = {}
for x in x_vals:
    x_val_avg_y[x] = 0
for i, x in enumerate(x_vals):
    x_val_avg_y[x] += y_vals[i]

for x in set(x_vals):
    x_val_avg_y[x] /= len([y for y in x_vals if y == x])

# Now approximating nutrient decay
nutrient_concentration_x = set(x_vals) # Gives us unique x values
nutrient_dict = {} # Initialize our dictionary

for x in nutrient_concentration_x:
    if x > 40:
        nutrient_dict[x] = 0.001 # Here, we say that if x > 40, then our ro is equal to 0.001. We set it 
                                 # to a non-zero value to prevent log issues later on. The cutoff is 40
                                 # as that is what we used earlier for calculating mu. It's when ro is very small.
                                 # It's a little arbitrary, but we're just taking guesses.
    else:
        nutrient_dict[x] = 0.2 - float(x_val_avg_y[x]) * a_guess # Otherwise, we approximate ro for each x value
                                                                 # by looking at the average y value for that x value
                                                                 # and using our "a_guess" to approximate how much nutrient
                                                                 # had been used to create the bacteria up to that point.

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
y_vals_log = np.array([pop_log_np, row_log_np]).transpose()

x_np = np.array(list(x_val_avg_y.keys()))

g_max, k, mu, a = curve_fit(log_diff, x_np, y_vals_log.flatten(), p0=(g_max_guess, k_guess, m1, a_guess))[0]
print(g_max, k, mu, a)
print(g_max_guess, k_guess, m1, a_guess)

# ODEINT extrapolation to 150 hours using fitted parameters
t = np.arange(4.5,150.1,0.1)
def modeling(x, t):
    n, row = x
    dndp = [(n*g_max*row)/(row+k) - (n*mu), (-a)*(n*g_max*row)/(row+k)]
    return dndp

result_1 = odeint(modeling, init_y, t)

# Fitting a 5th order polynomial
A = np.vstack([np.ones(len(x_vals)), x_vals, x_vals**2, x_vals**3, x_vals**4, x_vals**5]).T
print(A.shape, len(y_vals))

a, b, c, d, e, f = np.linalg.lstsq(A, y_vals, rcond=None)[0] #least squares solution to the linear matrix eqn.
xplot = np.arange(-10, 110)

plt.figure()
plt.plot(x_vals, y_vals, 'ro', label='Original data', markersize=4)
plt.plot(xplot, a+b*xplot+c*xplot**2+d*xplot**3+e*xplot**4+f*xplot**5, 'b', label='Fitted line')
plt.legend(loc=1)
plt.show()


# Fitting a 10th order polynomial
A2 = np.vstack([np.ones(len(x_vals)), x_vals, x_vals**2, x_vals**3, x_vals**4, x_vals**5, x_vals**6, x_vals**7, x_vals**8, x_vals**9, x_vals**10]).T 

a2, b2, c2, d2, e2, f2, g2, h2, i2, j2, k2 = np.linalg.lstsq(A2, y_vals, rcond=None)[0] #least squares solution to linear matrix eqn.

plt.figure()
plt.plot(x_vals, y_vals, 'ro', label='Original data', markersize=4)
plt.plot(xplot, a2+b2*xplot+c2*xplot**2+d2*xplot**3+e2*xplot**4+f2*xplot**5+g2*xplot**6+h2*xplot**7+i2*xplot**8+j2*xplot**9+k2*xplot**10, 'b', label='Fitted line')
plt.legend(loc=1)
plt.show()




# Plotting Extrapolation Graphs
plt.figure()
plt.plot(t, result_1, label = "ODEINT")
plt.plot(x_vals, y_vals, 'ro', label='Original data', markersize=4)
plt.xlabel("Hours")
plt.ylabel("Bacteria Population (CFU/ml)")
plt.title("ODEINT Extrapolation to 150 Hours")

plt.figure()
plt.plot(t, a+b*t+c*t**2+d*t**3+e*t**4+f*t**5, 'b', label='Polynomial')
plt.plot(x_vals, y_vals, 'ro', label='Original data', markersize=4)
plt.xlabel("Hours")
plt.ylabel("Bacteria Population (CFU/ml)")
plt.title("Polynomial Extrapolation to 150 Hours")