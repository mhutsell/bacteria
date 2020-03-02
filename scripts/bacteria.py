#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 14:39:43 2020

@authors: Mack, Camilla, Matthew
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as odeint

y_val_string = """3.00E+01
4.00E+01
2.00E+01
3.00E+01
1.20E+02
1.25E+02
2.35E+02
2.45E+02
3.65E+02
3.30E+02
7.85E+02
8.20E+02
8.84E+03
3.88E+03
3.70E+03
1.55E+05
1.77E+04
1.90E+04
3.00E+04
8.70E+04
7.04E+04
3.10E+05
5.00E+05
5.00E+05
5.50E+07
4.20E+07
5.50E+07
4.26E+07
1.13E+08
7.58E+07
1.26E+08
1.53E+08
7.95E+07
9.95E+07
1.19E+08
5.32E+07
1.98E+08
8.02E+07
1.29E+08
1.17E+08
1.14E+08
9.70E+07
9.50E+07
7.90E+07
9.44E+07
8.95E+07
2.40E+02
2.90E+02
4.40E+02
5.90E+02
1.29E+03
6.80E+02"""

x_val_string = """4
4
5
5
6
6
7
7
8
8
9
9
11
11
11
13.03
13.03
13.03
14.56
14.56
14.56
18.22
18.22
18.22
24
24
24
22
22
24
24
25.5
25.5
35
35
38
38
41
41
60
64
64
78
78
100
100
6
6
7
7
8
8
"""


y_vals = [float(y) for y in y_val_string.split()]
x_vals = [float(x) for x in x_val_string.split()]
value_dict = {}
for i, x in enumerate(x_vals):
    value_dict[x] = y_vals[i] 
    
x_vals_sorted = sorted(x_vals)
y_vals_sorted = [value_dict[x] for x in x_vals_sorted]
y_vals_log = [np.log(y) for y in y_vals_sorted]
plt.plot(x_vals_sorted, y_vals_sorted, 'ro')
plt.figure()
plt.plot(x_vals_sorted, y_vals_log, 'ro')
