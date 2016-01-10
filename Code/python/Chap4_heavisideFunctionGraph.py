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

H1 = lambda x, t: 0.5 + 1 / np.pi * np.arctan(x / t)
H2 = lambda x, t: 1 / (1 + np.exp(-x / t))
H3 = lambda x, t: np.exp(-np.exp(-x / t))
H4 = lambda x, t: 0.5 * (1 + np.tanh(x / t))

x = np.linspace(-10, 10, 1000)

# fileName = 'heaviside_function_example.eps'
# plt.figure(figsize=(30, 15))
# # plt.figure()
# plt.plot(x, H1(x, 1),
#          x, H2(x, 0.5),
#          x, H3(x, 1),
#          x, H4(x, 4),
#          x, S(x),
#          lw=linewidth)
# plt.xlabel(r'$x$')
# plt.ylabel(r'$\mathcal{H}(x)$')
# plt.grid('on')
# plt.minorticks_on()
# plt.grid(which='minor')
# plt.legend([r'$\mathcal{H}_1$', r'$\mathcal{H}_2$', r'$\mathcal{H}_3$', r'$\mathcal{H}_4$', r'$\mathcal{S}$'],
#            loc='best')
# plt.savefig(fileName, format='eps', dpi=1000, bbox_inches='tight')
# plt.show()

fileName = 'heaviside_function_with_control.eps'
plt.figure(figsize=(30, 15))
# plt.figure()
plt.plot(x, H2(x, 5),
         x, H2(x, 1),
         x, H2(x, 0.2),
         x, H2(x, 0.04),
         x, S(x),
         lw=linewidth)
plt.xlabel(r'$x$')
plt.ylabel(r'$\mathcal{H}(x)$')
plt.grid('on')
plt.minorticks_on()
plt.grid(which='minor')
plt.legend([r'$\eta = 5$', r'$\eta = 1$', r'$\eta = 0.2$', r'$\eta = 0.04$', r'$\mathcal{S}$'],
           loc='best')
plt.savefig(fileName, format='eps', dpi=1000, bbox_inches='tight')
plt.show()