__author__ = 'koorosh'
import numpy as np
import matplotlib.pyplot as plt
import pdb
# -----------------------------
font = {'family' : 'monospace',
        'weight' : 'normal',
        'size'   : 22}
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

def D(x, eta=1.0, x0=0):
    y = np.zeros([len(x), 1])
    x = x - x0
    for ix in range(0, len(x)):
        y[ix] = (1 - (np.tanh(x[ix] / eta))**2.0) / (2.0 * eta)
    return y
def T(x, eta=1.0, x0=0):
    y = np.zeros([len(x), 1])
    x = x - x0
    for ix in range(0, len(x)):
        y[ix] = (np.tanh(x[ix] / eta)**2 - 1) * np.tanh(x[ix] / eta) / eta**2
    return y
# ----------------------------
nx = 60
x = np.linspace(0, 1, nx).reshape(-1, 1)
dx = x[2] - x[1]
dt = 1e-5

L = genL(nx) / dx**2

wallVelocity = 1
BC = np.zeros([nx - 2, 1])
BC[0] = wallVelocity

un = np.zeros([nx, 1])
un[0] = BC[0]

X0 = 0.726
Xn = X0 # Wall location

p = 0.99
eta0 = 4 * dx / (2 * np.arctanh(p))
delta = D(x[1:-1], eta=eta0, x0=Xn).reshape(-1, 1)
doublet = T(x[1:-1], eta=eta0, x0=Xn).reshape(-1, 1)
Ft = 0
f = 0
alpha = -100
beta = -100

# Sensitivity analysis
sBC = np.zeros([nx - 2, 1])
sFt = 0
sf = 0
sun = np.zeros([nx, 1])
sf1 = 0
sf2 = 0
sf5 = 0
# Governing equations
for it in range(1, 50000):
    RHS1 = BC * dt / (2 * dx**2)
    RHS2 = np.eye(nx - 2, nx - 2).dot(un[1:-1])
    RHS3 = L[1:-1, :].dot(un) * dt / 2
    RHS = RHS1 + RHS2 + RHS3 + f
    A = np.eye(nx - 2, nx - 2) - L[1:-1, 1:-1] * dt / 2
    unp1 = np.linalg.solve(A, RHS)

    Un = np.trapz(delta * un[1:-1], x[1:-1], axis=0)
    Unp1 = np.trapz(delta * unp1, x[1:-1], axis=0)
    Ft = Ft + (Un + Unp1) * dt / 2
    Fc = Un
    F = alpha * Ft + beta * Fc
    f = F * delta * dt

    # ====== Sensitivity Analysis ====== #
    # ============== BEGIN ============= #
    sRHS1 = sBC * dt / (2 * dx**2)
    sRHS2 = np.eye(nx - 2, nx - 2).dot(sun[1:-1])
    sRHS3 = L[1:-1, :].dot(sun) * dt / 2
    sRHS = sRHS1 + sRHS2 + sRHS3 + sf
    A = np.eye(nx - 2, nx - 2) - L[1:-1, 1:-1] * dt / 2
    sunp1 = np.linalg.solve(A, sRHS)

    sf1t1 = np.trapz(delta * sun[1:-1], x[1:-1], axis=0)
    sf1t2 = np.trapz(delta * sunp1, x[1:-1], axis=0)
    sf1 = sf1 + (sf1t1 + sf1t2) * dt / 2

    sf2t1 = -np.trapz(doublet * un[1:-1], x[1:-1], axis=0)
    sf2t2 = -np.trapz(doublet * unp1, x[1:-1], axis=0)
    sf2 = sf2 + (sf2t1 + sf2t2) * dt / 2

    sf3 = np.trapz(delta * sun[1:-1], x[1:-1], axis=0)
    sf4 = -np.trapz(doublet * un[1:-1], x[1:-1], axis=0)

    sf5t1 = np.trapz(delta * un[1:-1], x[1:-1], axis=0)
    sf5t2 = np.trapz(delta * unp1, x[1:-1], axis=0)
    sf5 = sf5 + (sf5t1 + sf5t2) * dt / 2

    sf6 = np.trapz(delta * un[1:-1], x[1:-1], axis=0)

    sf = (alpha * (sf1 + sf2) + beta * (sf3 + sf4)) * delta - \
         (alpha * sf5 + beta * sf6) * doublet
    sf = sf * dt
    # ============== END =============== #
    # ====== Sensitivity Analysis ====== #
    un = np.append(np.append(BC[0], unp1), 0).reshape(-1, 1)
    sun = np.append(np.append(sBC[0], sunp1), 0).reshape(-1, 1)

xInd = np.argmax(x > Xn) - 1
uAnal = -BC[0] * x / Xn + BC[0]
suAnal = wallVelocity * x / Xn**2

rmsd = np.linalg.norm(un[:xInd] - uAnal[:xInd]) / (np.sqrt(len(un[:xInd])) * (np.max(np.abs(uAnal[:xInd])) -
                                                                              np.min(np.abs(uAnal[:xInd]))))
print(rmsd)

rmsd = np.linalg.norm(sun[:xInd] - suAnal[:xInd]) / (np.sqrt(len(sun[:xInd])) * (np.max(np.abs(suAnal[:xInd])) -
                                                                                 np.min(np.abs(suAnal[:xInd]))))
print(rmsd)

plt.figure()
plt.plot(x, sun,
         x, wallVelocity * x / X0**2)
plt.xlim([0, X0])
plt.figure()
plt.plot(x, un)
plt.show()
