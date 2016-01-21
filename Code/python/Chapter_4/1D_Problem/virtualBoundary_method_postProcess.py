__author__ = 'koorosh'
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sop
# -----------------------------
# plt.rc('text', usetex=True)
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

# ------------------------
n = np.array([5, 10, 15, 20, 40, 80])
RSME_RD = np.array([0.0413822258, 0.018991141, 0.0086346858, 0.0006043935, 0.0001588365, 0.003162531])
RSME_D = np.array([0.027840416, 0.0193350061, 0.0024083049, 0.0013924861, 0.0010884586, 0.0031240485])

fileName = 'effect_of_RD_on_simulation_vs_numberOfNodes_1D_problem.eps'
plt.figure(figsize=(30, 15))
plt.semilogy(n, RSME_RD, 'k',
             n, RSME_D, 'r--',
             lw=linewidth)
plt.legend(['RD function', 'Delta function'])
plt.xlabel('Number of nodes')
plt.ylabel('RSME')
plt.grid('on')
plt.minorticks_on()
plt.grid(which='minor')
plt.savefig(fileName, format='eps', dpi=1000, bbox_inches='tight')
plt.show()
# ------------------------