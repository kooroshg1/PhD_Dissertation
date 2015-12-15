__author__ = 'koorosh'
import numpy as np
import matplotlib.pyplot as plt

def phi1(r):
    y = np.zeros(len(r))
    for ir in range(0, len(r)):
        if np.abs(r[ir]) <= 1:
            y[ir] = 1 - np.abs(r[ir])
    return y

def phi2(r):
    y = np.zeros(len(r))
    for ir in range(0, len(r)):
        if np.abs(r[ir]) <= 0.5:
            y[ir] = (1 + np.sqrt(-3 * r[ir]**2 + 1)) / 3
        elif np.abs(r[ir]) <= 1.5:
            y[ir] = (5 - 3 * np.abs(r[ir]) - np.sqrt(-3 * (1 - np.abs(r[ir]))**2 + 1)) / 6
    return y

def phi3(r):
    y = np.zeros(len(r))
    for ir in range(0, len(r)):
        if np.abs(r[ir]) <= 1.0:
            y[ir] = (3 - 2 * np.abs(r[ir]) + np.sqrt(1 + 4 * np.abs(r[ir]) - 4 * r[ir]**2)) / 8
        elif np.abs(r[ir]) <= 1.5:
            y[ir] = (5 - 2 * np.abs(r[ir]) - np.sqrt(-7 + 12 * np.abs(r[ir]) - 4 * r[ir]**2)) / 8
    return y

def phi4(r):
    y = np.zeros(len(r))
    for ir in range(0, len(r)):
        if np.abs(r[ir]) <= 1.0:
            y[ir] = 61.0 / 112.0 \
                    - 11/42 * np.abs(r[ir]) \
                    - 11/56 * np.abs(r[ir])**2 \
                    + 1/12 * np.abs(r[ir])**3 \
                    + np.sqrt(3)/336 * np.sqrt(243 + 1584 * np.abs(r[ir]) - 784 * np.abs(r[ir])**2
                                               - 1560 * np.abs(r[ir])**3 + 500 * np.abs(r[ir])**4
                                               + 336 * np.abs(r[ir])**5 - 112 * np.abs(r[ir])**6)
        elif np.abs(r[ir]) <= 2.0:
            y[ir] = 21/16 + 7/12 * np.abs(r[ir]) - 7/8 * np.abs(r[ir])**2 + 1/6 * np.abs(r[ir])**3 - \
                    3/2 * phi4(np.abs(r[ir]) - 1)
        elif np.abs(r[ir]) <= 3:
            y[ir] = 9/8 - 23/12 * np.abs(r[ir]) + 3/4 * np.abs(r[ir])**2 - 1/12 * np.abs(r[ir])**3 + \
                    1/2 * phi4(np.abs(r[ir]) - 2)
    return y


def delta(x, x0, type):
    dx = x[2] - x[1]
    r = (x - x0) / dx
    out = phi4(r) / dx
    return out

x = np.linspace(-3, 3, 100)
# print(np.trapz(delta(x, 1, 2), x))
plt.figure()
plt.plot(x, delta(x, 1.3541, 2))
plt.show()
