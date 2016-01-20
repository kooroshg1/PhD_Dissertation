__author__ = 'koorosh'
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sop

# -----------------------------
font = {'family' : 'monospace',
        'weight' : 'normal',
        'size'   : 52}
plt.rc('font', **font)
linewidth = 5.0
markersize = 15
# font = {'family' : 'monospace',
#         'weight' : 'normal',
#         'size'   : 22}
# plt.rc('font', **font)
# linewidth = 2.0
# markersize = 15
# -----------------------------

# -----------------------------
# # RSME for RH vs. number of nodes
# n = np.array([11, 41, 81, 161])
# RSME_RH = np.array([0.01571, 0.0097, 0.0076, 0.0082])
# RSME_H = np.array([0.0474, 0.0177, 0.0077, 0.0092])
#
# fileName = 'effect_of_RH_on_simulation_vs_numberOfNodes_1D_problem.eps'
# plt.figure(figsize=(30, 15))
# plt.semilogy(n, RSME_RH, 'k',
#              n, RSME_H, 'r--',
#              lw=linewidth)
# plt.legend(['RH function', 'Step function'])
# plt.xlabel('Number of nodes')
# plt.ylabel('RSME')
# plt.grid('on')
# plt.minorticks_on()
# plt.grid(which='minor')
# plt.savefig(fileName, format='eps', dpi=1000, bbox_inches='tight')
# plt.show()
# -----------------------------

# -----------------------------
# RSME for RH vs. number of nodes
n = np.array([20, 40, 80, 120, 160, 180, 320])

# # x_wall = 0.4321, U = 1
# RSME_GE = np.array([0.1071337784, 0.0330065688, 0.0060771245, 0.004709103, 0.0030466798, 0.0038251288, 0.0046022806])
# RSME_SA = np.array([0.1371056254, 0.0356610516, 0.0183384031, 0.0063932112, 0.0046022806, 0.0057421065, 0.0046022806])

# # x_wall = 0.7583, U = 1
RSME_GE = np.array([0.0847649007, 0.0256031798, 0.0046463731, 0.003475902, 0.0023105673, 0.0028755434, 0.0044625426])
RSME_SA = np.array([0.1713840201, 0.0539913891, 0.0102838446, 0.0107751846, 0.0053758169, 0.0054250107, 0.0070518279])



def func(x, a, b):
        return a * x**b
popt, pcov = sop.curve_fit(f=func, xdata=n, ydata=RSME_GE)
aGE = popt[0]
bGE = popt[1]
print(bGE)
popt, pcov = sop.curve_fit(f=func, xdata=n, ydata=RSME_SA)
aSA = popt[0]
bSA = popt[1]
print(bSA)

fileName = 'penalizationMethod_SA_1D_problem_xw07583.eps'
plt.figure(figsize=(30, 15))
plt.loglog(n, RSME_GE, 'k',
           n, RSME_SA, 'r--',
           lw=linewidth)
plt.legend(['GE', 'SA'])
plt.xlabel('Number of nodes')
plt.ylabel('RSME')
plt.grid('on')
plt.minorticks_on()
plt.grid(which='minor')
plt.hold('on')
plt.loglog(n, aGE*n**bGE, 'k',
           n, aSA*n**bSA, 'r--', lw=1)
plt.savefig(fileName, format='eps', dpi=1000, bbox_inches='tight')
plt.show()

# -----------------------------

