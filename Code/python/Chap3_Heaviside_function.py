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
def H(x, x0=0):
    y = np.zeros([len(x), 1])
    for ix in range(0, len(x)):
        if x[ix] >= 0:
            y[ix] = 1
    return y

x = np.linspace(-2, 2, 1000)
y = H(x, x0=0)

plt.figure(figsize=(30, 15))
plt.plot(x, y, 'k',
         lw=linewidth)
plt.grid('on')
plt.xlabel(r'$x$')
plt.ylabel(r'$\mathcal{H}(x)$')
plt.savefig('Heaviside_Function7.eps', format='eps', dpi=1000, bbox_inches='tight')
plt.show()