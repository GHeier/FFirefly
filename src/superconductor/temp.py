import numpy as np
import matplotlib.pyplot as plt

wc = 0.1

def square_well_approximate(w):
    delta = 0.01
    return - 0.5 * (np.tanh((w + wc)/delta) - np.tanh((w - wc)/delta))

plt.figure()
warr = np.linspace(-0.5, 0.5, 1000)
farr = square_well_approximate(warr)
plt.plot(warr, farr)
plt.show()

