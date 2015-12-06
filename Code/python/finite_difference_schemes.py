__author__ = 'koorosh'
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
n = 8
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
F[-1] = Tn

T = np.linalg.solve(L, -F)
T = np.append(T, Tn)
T = np.append(T0, T)

# plt.figure(figsize=(30,15))
# plt.xlabel('x')
# plt.ylabel('T')
# plt.plot(x, Tc, 'k-',
#          x, x, 'wo',
#          ms=20, mew=7, mec='r', lw=linewidth)
# plt.legend(['FD', 'Analytical'], loc='best')
# plt.savefig('finitedifference_vs_analytical.eps', format='eps', dpi=1000, bbox_inches='tight')
# plt.show()

# # Sensitivity analysis

dDeltadL = 1 / (1.5 *n - 4)
Delta = x[2] - x[1]

T = T[1:-1]
dLdL = np.zeros([n-2, n-2])
dLdL[-1, -2] = -3
dLdL[-1, -1] = 1
dLdL = 0.5 * dDeltadL * 1 / Delta * dLdL
dFdL = 0.5 * dDeltadL * 1 / Delta * np.concatenate((np.zeros(n-3), [Tn]))
dFdL = dFdL - dLdL.dot(T)

dTdL = np.linalg.solve(L, dFdL)

print(T)
print(dTdL)
print(-x[1:-1]/(x[-1] - x[0])**2)

# print(np.divide(dTdL, (-x[1:-1]/(x[-1] - x[0])**2)))