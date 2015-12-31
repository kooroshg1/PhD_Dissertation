__author__ = 'koorosh'
import numpy as np
import matplotlib.pyplot as plt
# -----------------------------
font = {'family' : 'monospace',
        'weight' : 'normal',
        'size'   : 52}
plt.rc('font', **font)
linewidth = 9.0
markersize = 20
# -----------------------------
n = 161
x = np.linspace(0, 1, n)
dx = x[2] - x[1]
T0 = 0
Tn = 1

# # 2nd Order Central difference
# Assembling stiffness matrix
L = np.zeros([n-2, n-2])
L[0, 0] = -2
L[0, 1] = 1
L[-1, -2] = 1
L[-1, -1] = -2
for iRow in range(1, len(L)-1):
    L[iRow, iRow-1] = 1
    L[iRow, iRow] = -2
    L[iRow, iRow+1] = 1
# Assembling boundary condition
F = np.zeros([n-2, 1])
F[0] = T0
F[-1] = -(Tn - T0) / (x[-1] - x[0])

dTdb = np.linalg.solve(L, -F)
dTdb = np.append(dTdb, -1)
dTdb = np.append(0, dTdb)

dTdL_DSA = dTdb
dTdL_Anal = -x
rmsd = np.sqrt(np.sum((dTdL_DSA - dTdL_Anal)**2) / len(dTdL_Anal)) / (np.max(np.abs(dTdL_Anal)) -
                                                                      np.min(np.abs(dTdL_Anal)))
print(rmsd)
skip = 8
fileName = 'CSA_n' + np.str(n) + '.eps'
plt.figure(figsize=(30,15))
plt.plot(x, -x, 'k',
         x[::skip], dTdb[::skip], 'wo',
         ms=20, mew=7, mec='r', lw=linewidth)
plt.legend(['CSA', 'Analytical'], loc='best')
plt.xlabel('X')
plt.ylabel('dT/db')
# plt.savefig('CSA_vs_analytical.eps', format='eps', dpi=1000, bbox_inches='tight')
plt.show()