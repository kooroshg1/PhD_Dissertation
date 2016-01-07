__author__ = 'koorosh'
import numpy as np
import matplotlib.pyplot as plt
# -----------------------------
font = {'family' : 'monospace',
        'weight' : 'normal',
        'size'   : 52}
plt.rc('font', **font)
linewidth = 5.0
markersize = 15
# -----------------------------
nEl = 4
L = 1.0
EA = 1.0
P = 1 / np.pi
kSpring = 10.0
def FEA(nEl = 32, L = 1, EA = 1.0, P = 1 / np.pi, kSpring = 10.0):
    nx = nEl + 1
    x = np.linspace(0, L, nx)
    T = lambda x: np.sin(np.pi * x / L)
    Tdot = lambda x: -np.pi * x / L**2 * np.cos(np.pi * x / L)

    l = x[2] - x[1]
    k = EA / l * np.array([[1, -1], [-1, 1]])
    kdot = -EA / l**2 * np.array([[1, -1], [-1, 1]])

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
    Udot = np.linalg.solve(K, Fdot - Kdot.dot(U))
    return[U, Udot]

def ANAL(nEl = 32, L = 1, EA = 1.0, P = 1 / np.pi, kSpring = 10.0):
    x = np.linspace(0, L, nEl + 1)
    Uanal = (1/np.pi**2 * np.sin(np.pi * x) + 1/(10 * np.pi)).reshape(-1, 1)
    UdotAnal = ((-np.pi * x * np.cos(np.pi * x) - np.pi * (5 * x - 6) / 5 + 2 * np.sin(np.pi * x)) /
                (np.pi**2)).reshape(-1, 1)
    return [Uanal, UdotAnal]

def calcNRMSD(Y, y):
    # Y: Analytical, y: numerical
    return np.sqrt(np.sum((y - Y)**2) / len(y)) / (np.max(Y) - np.min(Y))

nEl = 4 * 2**np.arange(0, 6, 1)
NRMSD_GE = np.zeros([len(nEl), 1])
NRMSD_SA = np.zeros([len(nEl), 1])

for i in range(0, len(nEl)):
    U, Udot = FEA(nEl[i])
    Uanal, UdotAnal = ANAL(nEl[i])

    NRMSD_GE[i] = calcNRMSD(Uanal, U)
    NRMSD_SA[i] = calcNRMSD(UdotAnal[:-1], Udot[:-1])

plt.figure()
plt.semilogy(nEl, NRMSD_GE, 'k',
             nEl, NRMSD_SA,'k--',
             lw=linewidth, mew=linewidth, ms=markersize)
plt.legend(['Governing equations', 'Sensitivity analysis'])
plt.grid('on')
plt.xlabel('Number of elements')
plt.ylabel('NRMSE')
plt.show()