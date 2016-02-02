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
markersize = 15
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

print(sym.latex(D1.subs({t:1})))
print(sym.latex(D2.subs({t:1})))
print(sym.latex(D3.subs({t:1})))
print(sym.latex(D4.subs({t:1})))

D1 = sym.lambdify([x, t], D1, "numpy")
D2 = sym.lambdify([x, t], D2, "numpy")
D3 = sym.lambdify([x, t], D3, "numpy")
D4 = sym.lambdify([x, t], D4, "numpy")

t1 = 1
t2 = 1
t3 = 1
t4 = 1
x = np.linspace(-10, 10, 1000)
xInt = np.linspace(-100, 100, 1000000)
# print(np.trapz(D1(xInt, t1), xInt))
# print(np.trapz(D2(xInt, t2), xInt))
# print(np.trapz(D3(xInt, t3), xInt))
# print(np.trapz(D4(xInt, t4), xInt))

# fileName = 'delta_function_example.eps'
# plt.figure(figsize=(30, 15))
# plt.plot(x, D1(x, t1),
#          x, D2(x, t2),
#          x, D3(x, t3),
#          x, D4(x, t4),
#          lw=linewidth)
# plt.legend([r'$\mathcal{D}_1$', r'$\mathcal{D}_2$', r'$\mathcal{D}_3$', r'$\mathcal{D}_4$'])
# plt.xlabel(r'$x$')
# plt.ylabel(r'$\mathcal{D}(x)$')
# plt.grid('on')
# plt.minorticks_on()
# plt.grid(which='minor')
# plt.savefig(fileName, format='eps', dpi=1000, bbox_inches='tight')
# plt.show()

# fileName = 'delta_function_with_control.eps'
# plt.figure(figsize=(30, 15))
# plt.plot(x, D4(x, 2),
#          x, D4(x, 1),
#          x, D4(x, 0.25),
#          x, D4(x, 0.1),
#          lw=linewidth)
# plt.legend([r'$\eta = 2$', r'$\eta = 1$', r'$\eta = 0.25$', r'$\eta = 0.1$'])
# plt.xlabel(r'$x$')
# plt.ylabel(r'$\mathcal{D}(x)$')
# plt.grid('on')
# plt.minorticks_on()
# plt.grid(which='minor')
# plt.savefig(fileName, format='eps', dpi=1000, bbox_inches='tight')
# plt.show()