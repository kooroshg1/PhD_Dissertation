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
# -----------------------------
nx = 5
x = np.linspace(0, 1, nx+1).reshape(-1, 1)
x = x[:-1]
dx = x[2] - x[1]
dt = 0.01

U0 = 1
Un = 0
u = np.zeros([nx, 1])
u[0] = U0
u[-1] = Un

L = genL(nx) / dx**2
# L = L[1:-1, :]
un = u
unp1 = np.zeros([nx, 1])
unp1[0] = U0
unp1[-1] = Un
err = 1.0

f = 0
X0 = 0.5
Xn = X0
delta = genDelta(x, Xn, type='2point').reshape(-1, 1)
alpha = -10
beta = -0
Ft = 0
counter = 0
Unm1 = 0

I = np.eye(nx, nx)
for it in range(1, 200):
    BC = np.zeros([nx - 2, 1])
    BC[0] = (1 - dt/2)
    RHS = np.dot(I[1:-1, :] + dt * 0.5 * L[1:-1, :], un) + f * dt + BC
    A = I[1:-1, 1:-1] - dt * 0.5 * L[1:-1, 1:-1]
    unp1 = np.linalg.solve(A, RHS)
    un = np.append(np.append(1, unp1), 0).reshape(-1,1)
    #
    # Un = np.trapz(delta * un, x, axis=0)
    # Unp1 = np.trapz(delta * unp1, x, axis=0)
    # Ft = Ft + dt * 0.5 * (Un + Unp1)
    # Fc = Un
    # F = alpha * Ft + beta * Fc
    # f = F * delta
    # un = unp1

#
plt.figure()
plt.plot(unp1)
plt.show()
    # Un = np.trapz(delta * un[1:-1], x[1:-1], axis=0)
    # Unp1 = np.trapz(delta * unp1[1:-1], x[1:-1], axis=0)
    # timeIntegration = timeIntegration + dt * 0.5 * (Un + Unp1)
    # pdb.set_trace()
# while err > 1e-8:
# for it in range(0, 1000000):
#     counter += 1
#     unp1[1:-1] = un[1:-1] + dt * (L.dot(un) + f)
#     Un = np.trapz(delta * un[1:-1], x[1:-1], axis=0)
#     timeIntegration = timeIntegration + 0.5 * (Un + Unm1) * dt
#     Unm1 = np.trapz(delta * un[1:-1], x[1:-1], axis=0)
#     F = alpha * timeIntegration + beta * (Un)
#     f = F * delta
#     # print(Un)
#     # if counter > 5000:
#     #     break
#     err = np.max(np.abs(unp1[1:-1] - un[1:-1]))
#     un[1:-1] = unp1[1:-1]
#
# skip = 3
# plt.figure(figsize=(30, 15))
# plt.plot(x, unp1, 'k',
#          x[::skip], -x[::skip]/X0+1, 'wo',
#          lw=linewidth, mew=linewidth, ms=markersize)
# plt.xlim([0, X0])
# plt.ylim([0, 1])
# plt.xlabel('X')
# plt.ylabel('Response (u)')
# plt.legend(['IB', 'Analytical'])
# plt.savefig('peskin_method_0817.eps', format='eps', dpi=1000, bbox_inches='tight')
#
# plt.figure(figsize=(30, 15))
# plt.plot(x, np.abs(np.divide(unp1+x/X0-1, -x/X0+1)), 'k',
#          lw=linewidth, mew=linewidth, ms=markersize)
# plt.xlim([0, X0])
# plt.xlabel('X')
# plt.ylabel('Percentage of Error')
# plt.savefig('err_peskin_method_0817.eps', format='eps', dpi=1000, bbox_inches='tight')
# plt.show()