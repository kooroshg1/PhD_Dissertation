__author__ = 'koorosh'
import sympy as sym

EA, k, L, P, x = sym.symbols('EA k L P x', real=True)

u = 1 / (EA * k * sym.pi**2) * \
    (sym.pi * (EA * (2 * L - P * sym.pi) + k * (L - P * sym.pi) * (L - x)) + k * L**2 * sym.sin(sym.pi * x / L))

# u = sym.simplify(u.subs({L:1, EA:1, k:10, P:1/sym.pi}))
u = sym.simplify(u.subs({EA:1, k:10, P:1/sym.pi}))
print(sym.pretty(u))
