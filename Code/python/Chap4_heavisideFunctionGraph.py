__author__ = 'koorosh'
import numpy as np
import matplotlib.pyplot as plt
# -----------------------------
font = {'family' : 'monospace',
        'weight' : 'normal',
        'size'   : 52}
plt.rc('font', **font)
linewidth = 5.0
markersize = 15
# -----------------------------
def S(x):
    y = np.zeros(len(x))
    for ix in range(0, len(x)):
        if x[ix] == 0:
            y[ix] = 0.5
        elif x[ix] > 0:
            y[ix] = 1
    return y

H = lambda x: 1 / (1 + np.exp(-x))
x = np.linspace(-10, 10, 1000)
fileName = 'heaviside_function_example.eps'

plt.figure(figsize=(30, 15))
plt.plot(x, H(x), 'k',
         x, S(x), 'r',
         lw=linewidth)
plt.xlabel(r'$x$')
plt.ylabel(r'$\mathcal{H}(x)$')
plt.grid('on')
plt.legend(['Heaviside func.', 'Step func.'], loc='best')
plt.savefig(fileName, format='eps', dpi=1000, bbox_inches='tight')
plt.show()
