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
    print(max(values))
    print(min(values))
    n = 1000
    print("Plotting 111 path")
    for i in range(1, n):
        q = np.array([1, 1, 1]) * i / n
        q = np.dot(BZ, q)
        val = field(q)
        x.append(q[0])
        y.append(val)
    plt.plot(x, y)
    plt.show()

plot_111()
