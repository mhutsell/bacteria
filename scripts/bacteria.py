#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 14:39:43 2020

@authors: Mack, Camilla, Matthew
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import random
import scipy.stats as sp
from scipy.integrate import odeint
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score

data = pd.read_excel("../data/Bacteria_data.xlsx") # The file path here will have to be the same on your computer or
                                                   # the code will have an error. I formatted my excel sheet to have
                                                   # two columns— one with x at the top and one with y at the top.
                                                   # And I sorted the data in the excel sheet so that I wouldn't have
                                                   # to do it here. The excel sheet will be attached in our email.
                                                   # The organization of the files will also be on github.

data.plot(kind = "scatter", x = "X", y="Y") # Scatter plot our experimental data.
plt.title("Given Data: Bacterial Concentration vs Time (Hours)")
plt.xlabel("Time (Hours)")
plt.ylabel("Bacterial Concentration (CFU/ml")

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

# Plot log of y values vs x values for data points 40-51.
plt.figure()
plt.plot(x_vals_1, y_vals_1, 'ro')
plt.title("Log Y Values vs X values 40-51")
plt.xlabel("Time (Hours)")
plt.ylabel("log(bacterial concentration) (CFU/ml)")

# Small note on graphs: Where calculating r^2 is desirable, I have written it as an option but not enabled it.
# it should only take removing a few '#'s and parentheses and quotes and such to get it working if you want.

# Plot guessed y values vs x values for x values 40-51. Guessed y values are guesses of log(bacterial concentration).
# Guesses are made using m1 and b1 as previously calculated.
plt.figure()
y_vals_guessed_at_1 = [m1*x + b1 for x in x_vals_1] # Guess Y values 40-51 based on our mu_guess.
plt.plot(x_vals_1, y_vals_guessed_at_1, 'ro', color = "blue")
#r2_1 = r2_score(y_vals_1, y_vals_guessed_at_1)
plt.title("Guessed Log Y Values vs X Values 40-51") #r^2 = {}".format(round(r2_1, 2)))
plt.xlabel("Time (Hours)")
plt.ylabel("log(bacterial concentration) (CFU/ml)")


m1 *= -1 # Our calculated m1 is negative, which is correct, but in our differential equation, we subtract
         # n(mu). Therefore, to avoid two negatives becoming a positive, I make m1 a positive value by multiplying
         # by negative 1. Mu units: (1/hour)
print("Mu guess: {}".format(m1))


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

# Here we plot two figures.
# The first is our actual log y values from x values 0-31.
# The second is our guessed log y values from x values 0-31
plt.figure()
plt.plot(x_vals_2, y_vals_2, 'ro')
plt.title("Log Y Values vs X Values 0-31")
plt.xlabel("Time (Hours)")
plt.ylabel("log(bacterial concentration) (CFU/ml)")

plt.figure()
y_vals_guessed_at_2 = [m2*x + b2 for x in x_vals_2]
plt.plot(x_vals_2, y_vals_guessed_at_2, 'ro', color = "blue")
#r2_2 = r2_score(y_vals_2, y_vals_guessed_at_2)
plt.title("Guessed Log Y Values vs X Values 0-31") #r^2 = {}".format(round(r2_2, 2)))
plt.xlabel("Time (Hours)")
plt.ylabel("log(bacterial concentration) (CFU/ml)")

g_max_guess = m2+m1 # m2 = gmax - mu, so to find our gmax guess, we must add our mu guess (m1)
print("Gmax Guess: {}".format(g_max_guess)) # Our gmax guess. Gmax units: 1/hour


#################### Section 3: Finding a
# We know that, disregarding death of bacteria as we're just approximating, our
# amount of nutrient, 0.2 microgram/ml, led to the creation of roughly 2e8
# bacteria. Therefore, it takes 0,2 microgram/ml to make 2e8 bacteria per ml.
# so we can approximate 'a' by taking 2e8 and dividing it by 0.2. This will
# tell us how many bacteria one microgram/ml creates.

a_guess = 0.2 / np.max(y_vals)
print("'a' Guess: {}".format(a_guess)) # Our 'a' guess. 'a' units: μg/mL


#################### Section 4: Approximating nutrient decay due to bacterial growth

######## Finding the average y value for each x value
# We do this because when construct our "remaining nutrient" dictionary,
# It made more sense to us to have an average y value for each x value,
# So that we didn't have wildly different amounts for remaining nutrients
# on the same x value. It also made using a dictionary possible, as individual
# keys were unique vs trying to have a dictionary with 2 '4' keys. 
x_val_avg_y = {}
for x in x_vals:
    x_val_avg_y[x] = 0 # Initialize all of our dictionary values to 0
for i, x in enumerate(x_vals):
    x_val_avg_y[x] += y_vals[i] # Add our y values (all the values for 4 are added together)
    
# At this point, the dictionary has only unique x values and each x value is connected to the sum of
# each of the y values. For example 4: 30 + 40.

for x in set(x_vals):
    x_val_avg_y[x] /= len([y for y in x_vals if y == x]) # We divide by how many times they measured data at this time
                                                         # For example, 4: (30+40)/2


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
                                                                 # Nutrient units: μg/mL

approx_conc_d = [(y_vals[i+1]- y_vals[i]) for i in range(0,51)] # approximate derivative of bact. concentration by subtracting the 
                                                                # current point from the following point.
                                                                # Really bad approximation as the distance between points isn't even,
                                                                # but it's just for the guess.
approx_conc_d +=[approx_conc_d[50]] # so that we can have 52 values
k_guess_list = []

for i in range(32,41):
    # We calculate a guess for k by taking our differential equation for bacterial concentration and moving k to be alone.
    # We do this for 9 points in the middle of our data where we expect that ro will not be too big or too small.
    k_guess_at_i = ((y_vals[i] * g_max_guess * nutrient_dict[x_vals[i]]) / (y_vals[i] * m1 + approx_conc_d[i]) - nutrient_dict[x_vals[i]]) 
    k_guess_list.append(k_guess_at_i)

# We then get our K guess by averaging all of the K guesses that we made within the 'for' loop above.
k_guess = np.average(k_guess_list) # K units: μg/mL
print("K Guess: {}".format(k_guess))

#################### Section 5: curve_fit

init_y = [x_val_avg_y[x_vals[0]], 0.2] # Our initial values for our data. ro = 0.2 (μg/mL),
                                       # bacterial concentration = 35 (CFU/ml)

# This is our differential equation for bacterial concentration and ro
def diff_n_p(x, t, gmax, k, mu, a):
    n, row = x # x is an array with bacterial concentration in the first position, and ro in the second position
    dndp = [(n*gmax*row)/(row+k) - (n*mu), (-a)*(n*gmax*row)/(row+k)]  # Our two differential equations.
                                                                        # The first is for bacterial concentration and the second
                                                                        # is for ro.
    return dndp

# We have this function so that we can take the log of whatever results we get.
def log_diff(t, g_max_guess, k_guess, m1, a_guess):
    # Use odeint and assign it to a variable called "result" before taking the log so that we can deal w/ values < 1.
    result = odeint(diff_n_p, init_y, t, args=(g_max_guess, k_guess, m1, a_guess)).flatten()
    for i, x in enumerate(result):
        if x <= 0:
            # Added this because we were getting errors when taking the np.log(0)
            result[i] = 0.001
    return np.log(result)

pop_log_np = np.log(np.array(list(x_val_avg_y.values()))) # log of average bacterial concentrations per unique x values
row_log_np = np.log(np.array(list(nutrient_dict.values()))) # log of guessed nutrient values
y_vals_log = np.array([pop_log_np, row_log_np]).transpose() # combine our two above lists into one array which we take
                                                            # transpose of so that we can put it into curve_fit.

# Array with all of our unique x vaues. 
x_np = np.array(list(x_val_avg_y.keys()))

# Call curve_fit, calling our log_diff function (which then calls odeint) and passing in our args which are each of our
# guesses. The function returns fitted values for each of our four parameters. 
g_max, k, mu, a = curve_fit(log_diff, x_np, y_vals_log.flatten(), p0=(g_max_guess, k_guess, m1, a_guess))[0]

# Print our fitted value below our guess values to compare
print("Guessed Values: g_max: {}, k: {}, mu: {}, a: {}".format(g_max_guess, k_guess, m1, a_guess))
print("Fitted Values: g_max: {}, k: {}, mu: {}, a: {}".format(g_max, k, mu, a))


################ Section 6 Extrapolation to 150 hours

# ODEINT extrapolation to 150 hours using fitted parameters
t = np.arange(4.5,150.1,0.1) # Our time list ranging from 4.5->150 hours.

# Our diff equation (same as the diff_n_p function above, except we don't have to pass in args) 
def modeling(x, t):
    n, row = x # same as in diff_n_p
    dndp = [(n*g_max*row)/(row+k) - (n*mu), (-a)*(n*g_max*row)/(row+k)] # same as in diff_n_p 
    return dndp

result_1 = odeint(modeling, init_y, t) # our resulting bact. conc. values and ro values given 
                                       # our initial y values as initialized earlier
bacterial_result = [x[0] for x in result_1] # so that we have just the bacterial concentration values

# Plotting our odeint resulting bacterial concentration results over original data points to make sure it's reasonable
plt.figure()
plt.plot(t[:956], bacterial_result[:956], label = "ODEINT") # so that we only plot until hours = 100
plt.plot(x_vals, y_vals, 'ro', label='Original data', markersize=4)

# corresponding_y_vals = {}
# cor_y_v = []
# for x in set(x_vals):    
#     corresponding_y_vals[x] = bacterial_result[int(x*10-40)]
#     y = len([f for f in x_vals if f == x])
#     cor_y_v += [corresponding_y_vals[x]] * y                         
# r2_3 = r2_score(y_vals, cor_y_v)

plt.xlabel("Time (Hours)")
plt.ylabel("Bacteria Population (CFU/ml)")
plt.title("ODEINT Model over Original Data")

# new_x_val = []
# for i, x in enumerate(x_vals):
#     new_x_val.append(x)
#     if x in list(x_vals[0:i]):
#         new_x_val[i] = x + 0.001 * random.random()
# new_x_vals = np.asarray(new_x_val)

# Fitting a 5th order polynomial
# creating our 52 x 6 matrix to apply least square regression to.
A = np.vstack([np.ones(len(x_vals)), x_vals, x_vals**2, x_vals**3, x_vals**4, x_vals**5]).T

g, b, c, d, e, f = np.linalg.lstsq(A, y_vals, rcond=None)[0] #least squares solution to the linear matrix eqn.
# We use g here because a has already been assigned to our fitted value for a. 

# We allow for some room around where our data points lie so that we can see part of the polynomial's behavior 
# outside of the data points. That's why we have x values from -10 to 110 
xplot = np.arange(3.5, 101,0.1)

# In this figure we plot our 5th-order polynomial fit over our original data points so that we can see how well it fit.
plt.figure()
plt.plot(x_vals, y_vals, 'ro', label='Original data', markersize=4)
plt.plot(xplot, g+b*xplot+c*xplot**2+d*xplot**3+e*xplot**4+f*xplot**5, 'b', label='Fitted line')
plt.legend(loc=1)
#y_val_poly_5 = [g+b*t+c*t**2+d*t**3+e*t**4+f*t**5 for t in x_vals]
#r2_4 = r2_score(y_vals, y_val_poly_5)
plt.xlabel("Time (Hours)")
plt.ylabel("Bacteria Population (CFU/ml)")
plt.title("5th-Order Polynomial Fit on Data") # r^2: {}".format(round(r2_4,3)))
plt.show()


# Fitting a 10th order polynomial. Process is the same as above so detailed annotations are omitted. 
A2 = np.vstack([np.ones(len(x_vals)), x_vals, x_vals**2, x_vals**3, x_vals**4, x_vals**5, x_vals**6, x_vals**7, x_vals**8, x_vals**9, x_vals**10]).T 

a2, b2, c2, d2, e2, f2, g2, h2, i2, j2, k2 = np.linalg.lstsq(A2, y_vals, rcond=None)[0] #least squares solution to linear matrix eqn.

# Plot our 10th-order polynomial over our original data points to see the fit.
plt.figure()
plt.plot(x_vals, y_vals, 'ro', label='Original data', markersize=4)
plt.plot(xplot, a2+b2*xplot+c2*xplot**2+d2*xplot**3+e2*xplot**4+f2*xplot**5+g2*xplot**6+h2*xplot**7+i2*xplot**8+j2*xplot**9+k2*xplot**10, 'b', label='Fitted line')
plt.legend(loc=1)
#y_val_poly_10 = [a2+b2*xplot+c2*xplot**2+d2*xplot**3+e2*xplot**4+f2*xplot**5+g2*xplot**6+h2*xplot**7+i2*xplot**8+j2*xplot**9+k2*xplot**10 for xplot in x_vals]
#r2_5 = r2_score(y_vals, y_val_poly_10)
plt.xlabel("Time (Hours)")
plt.ylabel("Bacteria Population (CFU/ml)")
plt.title("10th-Order Polynomial Fit on Data") # r^2: {}".format(round(r2_5,3)))
plt.show()

# We did 5th order and 10th order just to show that both did poorly outside of the data range.
# From here on out we just use the coefficients from the 5th order polynomial.

# Plotting Extrapolation Graphs
# This first figure is our ODEINT model (over 150 hours) over our original data points to see if it behaves reasonably
# beyond the scope of our given data domain.


plt.figure()
plt.plot(t, bacterial_result, label = "ODEINT")
plt.plot(x_vals, y_vals, 'ro', label='Original data', markersize=4)
plt.xlabel("Time (Hours)")
plt.ylabel("Bacteria Population (CFU/ml)")
plt.title("ODEINT Extrapolation to 150 Hours")

# This second figure is our 5th-order polynomial model (over 150 hours) over our original data points to see how it behaves
# beyond the scope of our given data domain.
plt.figure()
plt.plot(t, g+b*t+c*t**2+d*t**3+e*t**4+f*t**5, 'b', label='Polynomial')
plt.plot(x_vals, y_vals, 'ro', label='Original data', markersize=4)
plt.xlabel("Time (Hours)")
plt.ylabel("Bacteria Population (CFU/ml)")
plt.title("Polynomial Extrapolation to 150 Hours")