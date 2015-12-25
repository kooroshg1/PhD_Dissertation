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

nx = 41
x = np.linspace(0, 1, nx + 2)
dx = x[2] - x[1]
dt = 1e-1

L = genL(nx)
Umwall = 1
Xswall = 0.842
indWn = np.argmax(x > Xswall) - 1
indWn1 = indWn - 1
un1w = 0
unw = (Xswall - x[indWn]) / (Xswall - x[indWn1]) * un1w

un = np.zeros(nx + 2)
I = np.eye(nx, nx)
If = np.zeros([nx, nx])
If[indWn - 1, indWn - 1] = 1
f = 0
un[0] = Umwall
for it in range(1, 100000):
    unp1 = I.dot(un[1:-1]) + dt * (L.dot(un) - f)
    unp1 = np.append(np.append(un[0], unp1), 0)
    unp1[indWn] = (Xswall - x[indWn]) / (Xswall - x[indWn1]) * unp1[indWn1]
    f = (unp1[1:-1] - un[1:-1]) / dt + L.dot(un)
    f = If.dot(f)
    un = list(unp1)


uAnal = - Umwall / Xswall * x + Umwall
uIB = un

xInd = np.argmax(x > Xswall) - 1
rmsd = np.sqrt(np.sum((uIB[:xInd] - uAnal[:xInd])**2) / len(uAnal)) / (np.max(np.abs(uAnal)) - np.min(np.abs(uAnal)))
print('RMSE = ', rmsd)

fileName = 'indirectForcing_wallLocation_' + np.str(np.int(Xswall*1000)) + '.eps'
skip = 1
plt.figure(figsize=(30, 15))
plt.plot(x, un, 'k',
         x[::skip], uAnal[::skip], 'wo',
         lw=linewidth, mew=linewidth, ms=markersize)
plt.xlim([0, Xswall])
plt.ylim([0, Umwall])
plt.xlabel('X')
plt.ylabel('Response (u)')
plt.legend(['IB', 'Analytical'])
plt.grid('on')
plt.savefig(fileName, format='eps', dpi=1000, bbox_inches='tight')
plt.show()