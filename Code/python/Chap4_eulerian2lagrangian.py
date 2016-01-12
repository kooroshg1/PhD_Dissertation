__author__ = 'koorosh'
import numpy as np
import matplotlib.pyplot as plt
import sympy as sym

# -----------------------------
font = {'family' : 'monospace',
        'weight' : 'normal',
        'size'   : 52}
plt.rc('font', **font)
linewidth = 5.0
markersize = 25
# -----------------------------
x, t = sym.symbols('x t', real=True)

H1 = 0.5 + 1 / sym.pi * sym.atan(x / t)
H2 = 1 / (1 + sym.exp(-x / t))
H3 = sym.exp(-sym.exp(-x / t))
H4 = 0.5 * (1 + sym.tanh(x / t))

D1 = sym.diff(H1, x)
D2 = sym.diff(H2, x)
D3 = sym.diff(H3, x)
D4 = sym.diff(H4, x)

D1 = sym.lambdify([x, t], D1, "numpy")
D2 = sym.lambdify([x, t], D2, "numpy")
D3 = sym.lambdify([x, t], D3, "numpy")
D4 = sym.lambdify([x, t], D4, "numpy")

x = np.linspace(0, 1, 9)
print(x)
y = lambda x: x**2
# X = 0.3321
# t = 0.1

X = 0.5813
t = 0.1

xt = x - X
Dbar1 = D1(x - X, t) / np.max(np.abs(D1(x - X, t)))
Dbar4 = D4(x - X, t) / np.max(np.abs(D4(x - X, t)))
print(np.trapz(D4(x - X, t), x))
print(np.trapz(np.multiply(y(x), D4(x - X, t)), x))
print(y(X))

fileName = 'euler2lagrange_example_03321.eps'
dx = x[2] - x[1]
# plt.figure()
# plt.plot(x, D4(xt, t))
# plt.figure()
# plt.plot(x, Dbar4)
plt.figure(figsize=(30, 15))
plt.plot(x, Dbar4 * np.trapz(np.multiply(y(x), D4(x - X, t)), x), 'g--s',
         x, y(x), 'k-o',
         X, np.trapz(np.multiply(y(x), D4(x - X, t)), x), 'ro',
         ms=markersize, lw=linewidth)
plt.hold('on')
plt.plot(X, y(X), 'b+',
         mew=linewidth, ms=markersize)
plt.xlabel(r'$x$')
plt.ylabel(r'$y(x)$')
plt.grid('on')
plt.minorticks_on()
plt.grid(which='minor')
# plt.savefig(fileName, format='eps', dpi=1000, bbox_inches='tight')
plt.show()