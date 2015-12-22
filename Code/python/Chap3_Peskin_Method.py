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
def phi2(r):
    y = np.zeros(len(r))
    for ir in range(0, len(r)):
        if np.abs(r[ir]) <= 1:
            y[ir] = 1 - np.abs(r[ir])
    return y
def phi3(r):
    y = np.zeros(len(r))
    for ir in range(0, len(r)):
        if np.abs(r[ir]) <= 0.5:
            y[ir] = (1 + np.sqrt(-3 * r[ir]**2 + 1)) / 3
        elif np.abs(r[ir]) <= 1.5:
            y[ir] = (5 - 3 * np.abs(r[ir]) - np.sqrt(-3 * (1 - np.abs(r[ir]))**2 + 1)) / 6
    return y
def phi4(r):
    y = np.zeros(len(r))
    for ir in range(0, len(r)):
        if np.abs(r[ir]) <= 1.0:
            y[ir] = (3 - 2 * np.abs(r[ir]) + np.sqrt(1 + 4 * np.abs(r[ir]) - 4 * r[ir]**2)) / 8
        elif np.abs(r[ir]) <= 1.5:
            y[ir] = (5 - 2 * np.abs(r[ir]) - np.sqrt(-7 + 12 * np.abs(r[ir]) - 4 * r[ir]**2)) / 8
    return y
def phi6(r):
    if type(r) is np.float64:
        r = np.array([r])
    y = np.zeros(len(r))
    for ir in range(0, len(r)):
        if np.abs(r[ir]) <= 1.0:
            y[ir] = 61.0 / 112.0 \
                    - 11/42 * np.abs(r[ir]) \
                    - 11/56 * np.abs(r[ir])**2 \
                    + 1/12 * np.abs(r[ir])**3 \
                    + np.sqrt(3)/336 * np.sqrt(243 + 1584 * np.abs(r[ir]) - 784 * np.abs(r[ir])**2
                                               - 1560 * np.abs(r[ir])**3 + 500 * np.abs(r[ir])**4
                                               + 336 * np.abs(r[ir])**5 - 112 * np.abs(r[ir])**6)
        elif np.abs(r[ir]) <= 2.0:
            y[ir] = 21/16 + 7/12 * np.abs(r[ir]) - 7/8 * np.abs(r[ir])**2 + 1/6 * np.abs(r[ir])**3 - \
                    3/2 * phi6(np.abs(r[ir]) - 1)
        elif np.abs(r[ir]) <= 3:
            y[ir] = 9/8 - 23/12 * np.abs(r[ir]) + 3/4 * np.abs(r[ir])**2 - 1/12 * np.abs(r[ir])**3 + \
                    1/2 * phi6(np.abs(r[ir]) - 2)
    return y
def genDelta(x, x0, type='2point'):
    dx = x[2] - x[1]
    r = (x - x0) / dx
    if type == '2point':
        out = phi2(r) / dx
    elif type == '3point':
        out = phi3(r) / dx
    elif type == '4point':
        out = phi4(r) / dx
    elif type == '6point':
        out = phi6(r) / dx
    else:
        print('Not available')
        out = 0
    return out
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
# -----------------------------
# ---------------------
nx = 161
x = np.linspace(0, 1, nx)
dx = x[2] - x[1]
dt = 1e-5

L = genL(nx) / dx**2

wallU = 10
BC = np.zeros([nx - 2, 1])
BC[0] = wallU

un = np.zeros([nx, 1])
un[0] = BC[0]

K = 10
X0 = 0.726
Xn = X0
deltaType = '2point'
delta = genDelta(x[1:-1], Xn, type=deltaType).reshape(-1, 1)

for it in range(1, 60000):
    RHS1 = BC * dt / (2 * dx**2)
    RHS2 = np.eye(nx - 2, nx - 2).dot(un[1:-1])
    RHS3 = L[1:-1, :].dot(un) * dt / 2
    U = np.trapz(delta * un[1:-1], x[1:-1], axis=0)
    Xnp1 = Xn + dt * U
    F = K * (Xnp1 - X0)
    f = F * delta * dt
    Xn = Xnp1
    RHS = RHS1 + RHS2 + RHS3 - f
    A = np.eye(nx - 2, nx - 2) - L[1:-1, 1:-1] * dt / 2
    unp1 = np.linalg.solve(A, RHS)
    err = np.max(np.abs(unp1 - un[1:-1]))
    if np.mod(it, 2000) == 0:
        print('Error = ', err)
    if err < 1e-6:
        print('Number of iterations = ', it)
        break
    un = np.append(np.append(BC[0], unp1), 0).reshape(-1, 1)

xInd = np.argmax(x > X0) - 1
uIB = un
uAnal = -BC[0] * x / X0 + BC[0]
rmsd = np.sqrt(np.sum((uIB[:xInd, 0] - uAnal[:xInd])**2) / len(uAnal)) / (np.max(np.abs(uAnal)) - np.min(np.abs(uAnal)))

print('RMSE = ', rmsd)

skip = 1
# fileName = 'classicalIB_wallVelocity_' + deltaType + '_' + np.str(K) + '.eps'
fileName = 'classicalIB_wallStiffness_' + deltaType + '_' + '10' + '.eps'
# print(fileName)
plt.figure(figsize=(30, 15))
plt.plot(x, un, 'k',
         x[::skip], -BC[0] * x[::skip] / X0 + BC[0], 'wo',
         lw=linewidth, mew=linewidth, ms=markersize)
plt.xlabel('X')
plt.ylabel('Response (u)')
plt.legend(['IB', 'Analytical'])
plt.xlim([0, X0])
plt.ylim([0, BC[0]])
plt.grid('on')
plt.savefig(fileName, format='eps', dpi=1000, bbox_inches='tight')
plt.show()

# #
# # # # 2nd Order Central difference
# # # Assembling stiffness matrix
# # def genL(nx):
# #     L = np.zeros([nx, nx])
# #     L[0, 0] = -2
# #     L[0, 1] = 1
# #     L[-1, -2] = 1
# #     L[-1, -1] = -2
# #     for iRow in range(1, len(L)-1):
# #         L[iRow, iRow-1] = 1
# #         L[iRow, iRow] = -2
# #         L[iRow, iRow+1] = 1
# #     return L
# #
# # L = genL(nx)
# # L = L[1:-1, :]
# # un = u
# # unp1 = np.zeros([nx, 1])
# # unp1[0] = U0
# # unp1[-1] = Un
# # err = 1.0
# # while err > 1e-5:
# #     unp1[1:-1] = un[1:-1] + dt * L.dot(un)
# #     err = np.max(np.abs(unp1[1:-1] - un[1:-1]))
# #     un[1:-1] = unp1[1:-1]
# #
# # plt.figure()
# # plt.plot(x, unp1)
# # plt.show()