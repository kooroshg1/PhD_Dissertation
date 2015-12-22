__author__ = 'koorosh'
import numpy as np
import matplotlib.pyplot as plt
import pdb
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
# ----------------------------
nx = 81
x = np.linspace(0, 1, nx).reshape(-1, 1)
dx = x[2] - x[1]
dt = 1e-5

L = genL(nx) / dx**2

wallVelocity = 10000
BC = np.zeros([nx - 2, 1])
BC[0] = wallVelocity

un = np.zeros([nx, 1])
un[0] = BC[0]

X0 = 0.726
Xn = X0 # Wall location
delta = genDelta(x[1:-1], Xn, type='2point').reshape(-1, 1) # Delta function
Ft = 0
f = 0
alpha = -1000
beta = -10
for it in range(1, 60000):
    RHS1 = BC * dt / (2 * dx**2)
    RHS2 = np.eye(nx - 2, nx - 2).dot(un[1:-1])
    RHS3 = L[1:-1, :].dot(un) * dt / 2
    RHS = RHS1 + RHS2 + RHS3 + f
    A = np.eye(nx - 2, nx - 2) - L[1:-1, 1:-1] * dt / 2
    unp1 = np.linalg.solve(A, RHS)
    err = np.max(np.abs(unp1 - un[1:-1]))
    if np.mod(it, 2000) == 0:
        print('Error = ', err)
    if err < 1e-8:
        print(it)
        break
    Un = np.trapz(delta * un[1:-1], x[1:-1], axis=0)
    Unp1 = np.trapz(delta * unp1, x[1:-1], axis=0)
    Ft = Ft + (Un + Unp1) * dt / 2
    Fc = Un
    F = alpha * Ft + beta * Fc
    f = F * delta * dt
    un = np.append(np.append(BC[0], unp1), 0).reshape(-1, 1)

xInd = np.argmax(x > X0) - 1
print('x @wall = ', x[xInd, 0])
uIB = un
uAnal = -BC[0] * x / X0 + BC[0]
rmsd = np.sqrt(np.sum((uIB[:xInd, 0] - uAnal[:xInd, 0])**2) / len(uAnal)) / (np.max(np.abs(uAnal)) -
                                                                             np.min(np.abs(uAnal)))

print('RMSE = ', rmsd)
print('U @wall = ', np.trapz(delta * un[1:-1], x[1:-1], axis=0))

fileName = 'virtualBoundary_wallVelocity_' + np.str(wallVelocity) + '.eps'
skip = 3
plt.figure(figsize=(30, 15))
plt.plot(x, un, 'k',
         x[::skip], -BC[0] * x[::skip] / Xn + BC[0], 'wo',
         lw=linewidth, mew=linewidth, ms=markersize)
plt.xlim([0, Xn])
plt.ylim([0, BC[0]])
plt.xlabel('X')
plt.ylabel('Response (u)')
plt.legend(['IB', 'Analytical'])
plt.grid('on')
plt.savefig(fileName, format='eps', dpi=1000, bbox_inches='tight')
# plt.savefig('vb_x06_u100.eps', format='eps', dpi=1000, bbox_inches='tight')
# plt.figure(figsize=(30, 15))
# plt.plot(x, un + BC[0] * x / Xn - BC[0], 'k',
#          lw=linewidth, mew=linewidth, ms=markersize)
# plt.xlim([0, Xn])
# plt.xlabel('X')
# plt.ylabel('Percentage of Error')
# plt.grid('on')
# plt.savefig('vb_x0385_u1_err.eps', format='eps', dpi=1000, bbox_inches='tight')
# plt.savefig('vb_x06_u100_err.eps', format='eps', dpi=1000, bbox_inches='tight')
plt.show()


# ---------------------
# Working implicit solver
# nx = 5
# x = np.linspace(0, 1, nx)
# dx = x[2] - x[1]
# dt = 0.1
#
# L = genL(nx) / dx**2
#
# BC = np.zeros([nx - 2, 1])
# BC[0] = 1
#
# un = np.zeros([nx, 1])
# un[0] = 1
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
#     un = np.append(np.append(1, unp1), 0).reshape(-1, 1)
#     print(un)





