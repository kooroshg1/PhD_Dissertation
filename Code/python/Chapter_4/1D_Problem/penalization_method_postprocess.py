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