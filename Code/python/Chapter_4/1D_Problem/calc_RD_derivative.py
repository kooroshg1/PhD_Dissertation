__author__ = 'koorosh'
import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
# -----------------------------
# plt.rc('text', usetex=True)
font = {'family' : 'monospace',
        'weight' : 'normal',
        'size'   : 52}
plt.rc('font', **font)
linewidth = 5.0
markersize = 15
# font = {'family' : 'monospace',
#         'weight' : 'normal',
#         'size'   : 22}
# plt.rc('font', **font)
# linewidth = 2.0
# markersize = 15
# -----------------------------
x, t = sym.symbols('x t', real=True)

H = 0.5 * (1 + sym.tanh(x / t))
D = sym.diff(H, x)
T = sym.diff(D, x)
print(sym.pretty(D))
print(sym.pretty(T))

D = sym.lambdify([x, t], D, "numpy")
T = sym.lambdify([x, t], T, "numpy")

x = np.linspace(-1.0, 1.0, 1000)

fileName = 'doubleFunction.eps'
plt.figure(figsize=(30, 15))
plt.plot(x, T(x, 0.5),
         x, T(x, 0.2),
         x, T(x, 0.1),
         lw=linewidth)
plt.legend([r'$\eta = 0.5$', r'$\eta = 0.2$', r'$\eta = 0.1$'])
plt.xlabel(r'$x$')
plt.ylabel(r'$T(x)$')
plt.grid('on')
plt.minorticks_on()
plt.grid(which='minor')
plt.savefig(fileName, format='eps', dpi=1000, bbox_inches='tight')
plt.show()
