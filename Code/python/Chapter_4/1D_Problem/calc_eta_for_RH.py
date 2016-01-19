__author__ = 'koorosh'
import numpy as np
import matplotlib.pyplot as plt

H = lambda x, eta: (1 + np.tanh(x / eta)) / 2.0

x = np.linspace(-3.0, 3.0, 100)
dx = x[1] - x[0]
dx = 0.01
p = 0.99
eta = dx / (2 * np.arctanh(p))
print(eta)
plt.figure()
plt.plot(x, H(x, eta),
         x, 0*x, 'ro')
plt.show()
