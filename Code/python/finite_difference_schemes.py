__author__ = 'koorosh gobal'
import numpy as np
import matplotlib.pyplot as plt
import sympy as sym

# n = 6
# x = np.linspace(0, 1, 6)
# dx = x[2] - x[1]
# T0 = 0
# Tn = 1
#
# # # 2nd Order Central difference
# # Assembling stiffness matrix
# Lc = np.zeros([n-2, n-2])
# Lc[0, 0] = -2
# Lc[0, 1] = 1
# Lc[-1, -2] = 1
# Lc[-1, -1] = -2
# for iRow in range(1, len(Lc)-1):
#     Lc[iRow, iRow-1] = 1
#     Lc[iRow, iRow] = -2
#     Lc[iRow, iRow+1] = 1
# # print(Lc)
# # Assembling boundary condition
# Fc = np.zeros([n-2, 1])
# Fc[0] = T0
# Fc[-1] = Tn
#
# T = np.linalg.solve(Lc, -Fc)
# T = np.append(T, Tn)
# T = np.append(T0, T)
# Tc = T
#

Di, Dip1, Din1, Din2 = sym.symbols('Di Dip1 Din1 Din2', real=True)
M = sym.Matrix([[Di, Di**2/2, Di**3/6, Di**4/24],
                [(Di+Dip1), (Di+Dip1)**2/2, (Di+Dip1)**3/6, (Di+Dip1)**4/24],
                [(-Din1), (-Din1)**2/2, (-Din1)**3/6, (-Din1)**4/24],
                [(-Din1-Din2), (-Din1-Din2)**2/2, (-Din1-Din2)**3/6, (-Din1-Din2)**4/24]])
Minv = M.inv().expand()
print(sym.pretty(Minv[-1, :]))