__author__ = 'koorosh'
import numpy as np
import matplotlib.pyplot as plt

def genL(nx):
    # # 2nd Order Central difference
    # Assembling stiffness matrix
    # L = np.zeros([nx, nx])
    # L[0, 0] = -2
    # L[0, 1] = 1
    # L[-1, -2] = 1
    # L[-1, -1] = -2
    # for iRow in range(1, len(L)-1):
    #     L[iRow, iRow-1] = 1
    #     L[iRow, iRow] = -2
    #     L[iRow, iRow+1] = 1
    # return L
    L = np.zeros([nx, nx + 2])
    for iRow in range(0, nx):
        L[iRow, iRow] = 1
        L[iRow, iRow + 1] = -2
        L[iRow, iRow + 2] = 1
    return L

nx = 4
x = np.linspace(0, 1, nx + 2)
dx = x[2] - x[1]
dt = 0.1

L = genL(nx)
Umwall = 1
Xswall = 0.5
indPhi = np.argmax(x > Xswall)
xswall = (Xswall - x[indPhi - 1]) / dx
Lp = np.zeros(L.shape)
Lline = np.zeros(nx + 2)
Lline[indPhi] = 1
Lline[indPhi - 1] = (1 - xswall) / xswall
print(L)
print()
L[indPhi - 1 - 1, :] = L[indPhi - 1 - 1, :] - Lline
L[indPhi + 1 - 1, :] = L[indPhi + 1 - 1, :] - Lline
print(np.round(L))
# BC = np.zeros([nx - 2, 1])
# BC[0] = Umwall
#
# un = np.zeros([nx, 1])
# un[0] = Umwall
# for it in range(1, 600):
#     RHS1 = BC * dt / (2 * dx**2)
#     RHS2 = np.eye(nx - 2, nx - 2).dot(un[1:-1])
#     RHS3 = L[1:-1, :].dot(un) * dt / 2
#     RHS = RHS1 + RHS2 + RHS3
#     A = np.eye(nx - 2, nx - 2) - L[1:-1, 1:-1] * dt / 2
#     unp1 = np.linalg.solve(A, RHS)
#     err = np.max(np.abs(unp1 - un[1:-1]))
#     if err < 1e-5:
#         break
#     un = np.append(np.append(Umwall, unp1), 0).reshape(-1, 1)
#
# print(un)
# plt.figure()
# plt.plot(x, un)
# plt.show()