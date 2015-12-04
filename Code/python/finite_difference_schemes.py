__author__ = 'koorosh gobal'
import numpy as np
import matplotlib.pyplot as plt
# -----------------------
font = {'family' : 'monospace',
        'weight' : 'normal',
        'size'   : 52}
plt.rc('font', **font)
linewidth = 9.0
markersize = 20
# -------------------------
n = 6
x = np.linspace(0, 1, 6)
dx = x[2] - x[1]
T0 = 0
Tn = 1

# # 2nd Order Central difference
# Assembling stiffness matrix
Lc = np.zeros([n-2, n-2])
Lc[0, 0] = -2
Lc[0, 1] = 1
Lc[-1, -2] = 1
Lc[-1, -1] = -2
for iRow in range(1, len(Lc)-1):
    Lc[iRow, iRow-1] = 1
    Lc[iRow, iRow] = -2
    Lc[iRow, iRow+1] = 1
# print(Lc)
# Assembling boundary condition
Fc = np.zeros([n-2, 1])
Fc[0] = T0
Fc[-1] = Tn

T = np.linalg.solve(Lc, -Fc)
T = np.append(T, Tn)
T = np.append(T0, T)
Tc = T

plt.figure(figsize=(30,15))
plt.xlabel('x')
plt.ylabel('T')
plt.plot(x, Tc, 'k-',
         x, x, 'wo',
         ms=20, mew=7, mec='r', lw=linewidth)
plt.legend(['FD', 'Analytical'], loc='best')
plt.savefig('finitedifference_vs_analytical.eps', format='eps', dpi=1000, bbox_inches='tight')
plt.show()