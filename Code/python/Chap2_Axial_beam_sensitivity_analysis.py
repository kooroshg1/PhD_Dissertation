__author__ = 'koorosh'
import numpy as np
import matplotlib.pyplot as plt

nEl = 4
nx = nEl + 1
L = 1.01
EA = 1.0
P = 1 / np.pi
x = np.linspace(0, L, nx)
T = lambda x: np.sin(np.pi * x / L)
Tdot = lambda x: -np.pi / L**2 * np.cos(np.pi * x / L)
kSpring = 10.0

l = x[2] - x[1]
k = EA / l * np.array([[1, -1], [-1, 1]])
kdot = -EA / (1/(nx-1)) * np.array([[1, -1], [-1, 1]])

# Assembling the stiffness matrix
K = 0
F = 0
Kdot = 0
Fdot = 0
for iEl in range(1, nEl + 1):
    nS, nE = [iEl - 1, iEl]
    xS, xE = np.multiply(l, [nS, nE])

    fMat = np.zeros([nx, 1])
    fMat[nS, 0] = 1.0
    fMat[nE, 0] = 1.0
    xEl = np.linspace(xS, xE, 10)
    f = np.trapz(T(xEl), xEl) / 2
    F = F + np.multiply(f, fMat)
    fdot = np.trapz(Tdot(xEl), xEl) / 2
    Fdot = Fdot + np.multiply(fdot, fMat)

    Kmat = np.zeros([nx, nx])
    Kmat[nS, nS] = 1
    Kmat[nE, nE] = 1
    Kmat = Kmat[:, [nS, nE]]
    K = K + Kmat.dot(k).dot(Kmat.T)
    if iEl == nEl:
        Kdot = Kdot + Kmat.dot(kdot).dot(Kmat.T)

# # Sanity check for distributed load
# print(np.trapz(T(np.linspace(0, L, 1000)), np.linspace(0, L, 1000)))
# print(np.sum(F))

# Model spring at the boundary
Kspring = np.zeros([nx, nx])
Kspring[-1, -1] = -kSpring
K = K - Kspring

# Add tip load to the F vector
F[0] = F[0] - P

U = np.linalg.solve(K, F)
Uanal = (1/np.pi**2 * np.sin(np.pi * x) + 1/(10 * np.pi)).reshape(-1,1)
NRMSD = np.sqrt(np.sum((U - Uanal)**2) / len(U)) / (np.max(Uanal) - np.min(Uanal))
print(NRMSD)

Udot = np.linalg.solve(K - Kspring, Fdot - Kdot.dot(U))
print(Udot)

# np.savetxt('U.txt', U)

# plt.figure()
# plt.plot(x, U, 'k',
#          x, Uanal, 'r--')
# plt.legend(['FEA', 'EXACT'])
# plt.show()