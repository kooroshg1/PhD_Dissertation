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

# # ------------------------
# n = np.array([5, 10, 15, 20, 40, 80])
# RSME_RD = np.array([0.0413822258, 0.018991141, 0.0086346858, 0.0006043935, 0.0001588365, 0.003162531])
# RSME_D = np.array([0.027840416, 0.0193350061, 0.0024083049, 0.0013924861, 0.0010884586, 0.0031240485])
#
# fileName = 'effect_of_RD_on_simulation_vs_numberOfNodes_1D_problem.eps'
# plt.figure(figsize=(30, 15))
# plt.semilogy(n, RSME_RD, 'k',
#              n, RSME_D, 'r--',
#              lw=linewidth)
# plt.legend(['RD function', 'Delta function'])
# plt.xlabel('Number of nodes')
# plt.ylabel('RSME')
# plt.grid('on')
# plt.minorticks_on()
# plt.grid(which='minor')
# plt.savefig(fileName, format='eps', dpi=1000, bbox_inches='tight')
# plt.show()
# # ------------------------

# # ------------------------
# n = np.array([20, 40, 60, 80, 100, 120, 160])
#
# # X0 = 0.6981 | U = 1.0
# RSME_GE = np.array([0.0330216527, 0.0147602907, 0.0095046884, 0.0070085343, 0.0055504215, 0.0045943183, 0.0034321797])
# RSME_SA = np.array([0.0486864679, 0.0297240499, 0.0233176565, 0.0202478105, 0.0185352937, 0.0175085386, 0.0066284337])
#
# # # X0 = 0.4327 | U = 1000.0
# # RSME_GE = np.array([0.0510339145, 0.0238020131, 0.0152753218, 0.0112147687, 0.0089620359, 0.0074083038, 0.0055204187])
# # RSME_SA = np.array([0.0771418053, 0.0428210249, 0.0293896996, 0.026107942, 0.0177694389, 0.0169613073, 0.0121552498])
#
# def func(x, a, b):
#         return a * x**b
# popt, pcov = sop.curve_fit(f=func, xdata=n, ydata=RSME_GE)
# aGE = popt[0]
# bGE = popt[1]
# print(bGE)
# popt, pcov = sop.curve_fit(f=func, xdata=n, ydata=RSME_SA)
# aSA = popt[0]
# bSA = popt[1]
# print(bSA)
#
# fileName = 'virtualBoundaryMethod_SA_1D_problem_xw06981_U1.eps'
# plt.figure(figsize=(30, 15))
# plt.loglog(n, RSME_GE, 'k',
#            n, RSME_SA, 'r--',
#            lw=linewidth)
# plt.legend(['GE', 'SA'])
# plt.xlabel('Number of nodes')
# plt.ylabel('RMSE')
# plt.grid('on')
# plt.minorticks_on()
# plt.grid(which='minor')
# plt.hold('on')
# plt.loglog(n, aGE*n**bGE, 'k',
#            n, aSA*n**bSA, 'r--', lw=1)
# plt.savefig(fileName, format='eps', dpi=1000, bbox_inches='tight')
# plt.show()
# # ------------------------

# # ------------------------
# # U = 1
# # Xn = 0.5321
# # x = np.loadtxt('dataFiles/x_VirtualB_05321_n95.txt')
# # dUdL_CSA = np.loadtxt('dataFiles/dUdL_VirtualB_05321_n95.txt')
#
# # U = 10
# # Xn = 0.6845
# # x = np.loadtxt('dataFiles/x_VirtualB_06845_n51.txt')
# # dUdL_CSA = np.loadtxt('dataFiles/dUdL_VirtualB_06845_n51.txt')
#
# # U = 100
# # Xn = 0.2861
# # x = np.loadtxt('dataFiles/x_VirtualB_02861_n71.txt')
# # dUdL_CSA = np.loadtxt('dataFiles/dUdL_VirtualB_02861_n71.txt')
#
# U = 1000
# Xn = 0.7412
# x = np.loadtxt('dataFiles/x_VirtualB_07412_n86.txt')
# dUdL_CSA = np.loadtxt('dataFiles/dUdL_VirtualB_07412_n86.txt')
#
# dUdL_Anal = U * x / Xn**2
#
# fileName = 'virtualBoundary_sensitivityProfile_xw0' + str(int(Xn * 10000)) + '_U' + str(int(U)) + '.eps'
# plt.figure(figsize=(30, 15))
# plt.plot(x, dUdL_CSA, 'wo',
#          x, dUdL_Anal, 'k',
#          ms=markersize, lw=linewidth, mew=linewidth)
# plt.legend(['CSA', 'Analytical'], loc='best')
# plt.xlim([0, Xn])
# plt.xlabel(r'$X$')
# plt.ylabel(r'$\partial u/\partial L$')
# plt.grid('on')
# plt.minorticks_on()
# plt.grid(which='minor')
# plt.savefig(fileName, format='eps', dpi=1000, bbox_inches='tight')
# plt.show()
# # ------------------------

# # ------------------------
# U = 1
# Xn = 0.5321
# x = np.loadtxt('dataFiles/x_VirtualB_05321_n95.txt')
# dUdL_CSA = np.loadtxt('dataFiles/dUdL_VirtualB_05321_n95.txt')

# U = 10
# Xn = 0.6845
# x = np.loadtxt('dataFiles/x_VirtualB_06845_n51.txt')
# dUdL_CSA = np.loadtxt('dataFiles/dUdL_VirtualB_06845_n51.txt')

# U = 100
# Xn = 0.2861
# x = np.loadtxt('dataFiles/x_VirtualB_02861_n71.txt')
# dUdL_CSA = np.loadtxt('dataFiles/dUdL_VirtualB_02861_n71.txt')

U = 1000
Xn = 0.7412
x = np.loadtxt('dataFiles/x_VirtualB_07412_n86.txt')
dUdL_CSA = np.loadtxt('dataFiles/dUdL_VirtualB_07412_n86.txt')

dUdL_Anal = U * x / Xn**2

wIndex = np.where(x > Xn)[0][0]
m = (dUdL_CSA[wIndex - 4] - dUdL_CSA[wIndex - 5]) / (x[wIndex - 4] - x[wIndex - 5])

dUdL_CSA[wIndex - 3] = m * (x[wIndex - 3] - x[wIndex - 4]) + dUdL_CSA[wIndex - 4]
dUdL_CSA[wIndex - 2] = m * (x[wIndex - 2] - x[wIndex - 4]) + dUdL_CSA[wIndex - 4]
dUdL_CSA[wIndex - 1] = m * (x[wIndex - 1] - x[wIndex - 4]) + dUdL_CSA[wIndex - 4]
dUdL_CSA[wIndex] = m * (x[wIndex] - x[wIndex - 4]) + dUdL_CSA[wIndex - 4]

fileName = 'virtualBoundary_sensitivityProfile_xw0' + str(int(Xn * 10000)) + '_U' + str(int(U)) + '_reconstructed.eps'
plt.figure(figsize=(30, 15))
plt.plot(x, dUdL_CSA, 'wo',
         x, dUdL_Anal, 'k',
         ms=markersize, lw=linewidth, mew=linewidth)
plt.legend(['CSA', 'Analytical'], loc='best')
plt.xlim([0, Xn])
plt.xlabel(r'$X$')
plt.ylabel(r'$\partial u/\partial L$')
plt.grid('on')
plt.minorticks_on()
plt.grid(which='minor')
plt.savefig(fileName, format='eps', dpi=1000, bbox_inches='tight')
plt.show()
# # ------------------------