__author__ = 'koorosh'
import numpy as np
import matplotlib.pyplot as plt

D = lambda x, eta: (-np.tanh(x / eta)**2 + 1) / (2 * eta)

x = np.linspace(-3, 3, 1000)

p = 0.99
R = 0.01
eta = R / (np.arctanh(np.sqrt(p)))
print(eta)
# plt.figure()
# plt.plot(x, D(x, eta))
# plt.show()