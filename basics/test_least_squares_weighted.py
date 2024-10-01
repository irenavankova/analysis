#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def func(x, a, b):

    return a * x + b


xdata = np.linspace(0, 4, 50)
w = np.zeros(len(xdata)) + 0.000001

y = func(xdata, 1, 0)

i_bad = np.array([1, 5, 10, 25])

y_bad = func(xdata, 1, 2)

ydata = y

ydata[i_bad] = y_bad[i_bad]
w[i_bad] = 0.9


plt.plot(xdata, ydata, 'bo', label='data')

popt, pcov = curve_fit(func, xdata, ydata)
plt.plot(xdata, func(xdata, *popt), 'k-', label='fit: a=%5.3f, b=%5.3f' % tuple(popt))

popt, pcov = curve_fit(func, xdata, ydata, sigma=w)
plt.plot(xdata, func(xdata, *popt), 'r--', label='fit: a=%5.3f, b=%5.3f' % tuple(popt))

plt.show()