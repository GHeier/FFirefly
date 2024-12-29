import fcode
import fcode.config as cfg

import numpy as np
import matplotlib.pyplot as plt

BZ = np.array(cfg.brillouin_zone)

def plot_111():
    x, y = list(), list()
    field = fcode.ScalarField("chi_mesh_static.dat")
    points = field.points
    values = field.values
    n = 100
    print("Plotting 111 path")
    for i in range(0, n):
        q = np.array([0.5, 0.5, 0.5]) * i / n
        q = np.dot(BZ, q)
        val = field(q)
        x.append(q[0])
        y.append(val)
    plt.plot(x, y)
    plt.show()

plot_111()
