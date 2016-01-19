__author__ = 'koorosh'
import numpy as np
import matplotlib.pyplot as plt
import time
# -----------------------------
# font = {'family' : 'monospace',
#         'weight' : 'normal',
#         'size'   : 52}
# plt.rc('font', **font)
# linewidth = 5.0
# markersize = 15
font = {'family' : 'monospace',
        'weight' : 'normal',
        'size'   : 22}
plt.rc('font', **font)
linewidth = 2.0
markersize = 15
# -----------------------------
def H(x, eta=1, x0=0):
    y = np.zeros([len(x), 1])
    for ix in range(0, len(x)):
        if x[ix] - x0 >= 0:
            y[ix] = 1
    return y

# def H(x, eta=1.0, x0=0):
#     y = np.zeros([len(x), 1])
#     x = x - x0
#     for ix in range(0, len(x)):
#         y[ix] = (1 + np.tanh(x[ix] / eta)) / 2.0
#     return y

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

nx = 81
x = np.linspace(0, 1, nx).reshape(-1, 1)
dx = x[2] - x[1]
dt = 1e-4

L = genL(nx) / dx**2

U = 1000
BC = np.zeros([nx - 2, 1])
BC[0] = U

un = np.zeros([nx, 1])
un[0] = U
kappa = 2/dt
Xn = 0.4327
eta0 = dx / (2 * np.arctanh(0.99))

xInd = np.argmax(x > Xn) - 1
uAnal = -BC[0] * x / Xn + BC[0]

# plt.figure()
for it in range(1, 40000):
    RHS1 = BC * dt / (2 * dx**2)
    RHS2 = np.eye(nx - 2, nx - 2).dot(un[1:-1])
    RHS3 = L[1:-1, :].dot(un) * dt / 2
    f = -kappa * np.multiply(un[1:-1], H(x[1:-1], eta=eta0, x0=Xn)) * dt
    # print(f)
    RHS = RHS1 + RHS2 + RHS3 + f
    A = np.eye(nx - 2, nx - 2) - L[1:-1, 1:-1] * dt / 2
    unp1 = np.linalg.solve(A, RHS)
    err = np.max(np.abs(unp1 - un[1:-1]))
    if err < 1e-10:
        print(it)
        break
    un = np.append(np.append(BC[0], unp1), 0).reshape(-1, 1)
    # plt.plot(x, un,
    #          x, -BC[0] * x / Xn + BC[0], 'r')
    # plt.title(str(it))
    # plt.xlim([0, Xn])
    # plt.ylim([0, BC[0]])
    # plt.pause(1e-10)
    # plt.clf()
    # rmsd = np.linalg.norm(un[:xInd] - uAnal[:xInd]) / (np.sqrt(len(un[:xInd])) * (np.max(np.abs(uAnal)) -
    #                                                                            np.min(np.abs(uAnal))))
    # print(rmsd)
    # print(un)

rmsd = np.linalg.norm(un[:xInd] - uAnal[:xInd]) / (np.sqrt(len(un[:xInd])) * (np.max(np.abs(uAnal)) -
                                                                              np.min(np.abs(uAnal))))
# np.savetxt('uIB_Hf_n161.txt', un)
# np.savetxt('uAnal.txt', uAnal)
print(rmsd)
# skip = 1
# # plt.figure()
# plt.figure(figsize=(30, 15))
# plt.plot(x, un, 'k',
#          x[::skip], -BC[0] * x[::skip] / Xn + BC[0], 'wo',
#          lw=linewidth, mew=linewidth, ms=markersize)
# plt.xlabel('X')
# plt.ylabel('Response (u)')
# plt.legend(['IB', 'Analytical'])
# plt.xlim([0, Xn])
# plt.ylim([0, U])
# plt.grid('on')
# plt.show()
# #
