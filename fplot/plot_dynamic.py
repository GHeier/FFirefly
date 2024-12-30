import fcode
import fcode.config as cfg
from fcode import *

import numpy as np
import matplotlib.pyplot as plt

BZ = np.array(cfg.brillouin_zone)

def plot_111():
    field = fcode.ScalarField("chi_mesh_dynamic.dat", 4, True)
    n = 100
    wpts = 1
    print("Plotting 111 path")
    for w in range(0, wpts):
        x, y = list(), list()
        for i in range(1, n):
            q = np.array([0.5, 0.5, 0.5]) * i / n
            q = np.dot(BZ, q)
            q = Vec(q[0], q[1], q[2], 0.0032)
            val = field(q)
            x.append(q.x)
            y.append(val)
        print(y)
        plt.plot(x, y, label="w = " + str(w / wpts))
    plt.show()

plot_111()

