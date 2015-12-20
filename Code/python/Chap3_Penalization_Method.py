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
def H(x, x0=0):
    y = np.zeros([len(x), 1])
    for ix in range(0, len(x)):
        if x[ix] - x0 >= 0:
            y[ix] = 1
    return y
def genL(nx):
    # # 2nd Order Central difference
    # Assembling stiffness matrix
    L = np.zeros([nx, nx])
    L[0, 0] = -2
    L[0, 1] = 1
    L[-1, -2] = 1
    L[-1, -1] = -2
    for iRow in range(1, len(L)-1):
        L[iRow, iRow-1] = 1
        L[iRow, iRow] = -2
        L[iRow, iRow+1] = 1
    return L

nx = 50
x = np.linspace(0, 1, nx).reshape(-1, 1)
dx = x[2] - x[1]
dt = 0.0001

L = genL(nx) / dx**2

BC = np.zeros([nx - 2, 1])
BC[0] = 1

un = np.zeros([nx, 1])
un[0] = 1
kappa = 10001
Xn = 0.65
for it in range(1, 60000):
    RHS1 = BC * dt / (2 * dx**2)
    RHS2 = np.eye(nx - 2, nx - 2).dot(un[1:-1])
    RHS3 = L[1:-1, :].dot(un) * dt / 2
    f = -kappa * np.multiply(un[1:-1], H(x[1:-1], x0=Xn)) * dt
    # print(f)
    RHS = RHS1 + RHS2 + RHS3 + f
    A = np.eye(nx - 2, nx - 2) - L[1:-1, 1:-1] * dt / 2
    unp1 = np.linalg.solve(A, RHS)
    err = np.max(np.abs(unp1 - un[1:-1]))
    if err < 1e-10:
        print(it)
        break
    un = np.append(np.append(1, unp1), 0).reshape(-1, 1)
    # print(un)

plt.figure()
# plt.figure(figsize=(30, 15))
plt.plot(x, un, 'k',
         x, -1 * x / Xn + 1, 'wo',
         lw=linewidth, mew=linewidth, ms=markersize)
plt.xlabel('X')
plt.ylabel('Response (u)')
plt.legend(['IB', 'Analytical'])
plt.xlim([0, Xn])
plt.ylim([0, 1])
plt.grid('on')
# plt.savefig('vb_x0385_u1.eps', format='eps', dpi=1000, bbox_inches='tight')

plt.figure()
# plt.figure(figsize=(30, 15))
plt.plot(x, un + 1 * x / Xn - 1, 'k',
         lw=linewidth, mew=linewidth, ms=markersize)
plt.xlim([0, Xn])
plt.xlabel('X')
plt.ylabel('Percentage of Error')
plt.grid('on')
# plt.savefig('vb_x0385_u1_err.eps', format='eps', dpi=1000, bbox_inches='tight')
plt.show()


