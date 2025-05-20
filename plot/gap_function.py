import firefly
import firefly.config as cfg

import numpy as np
import matplotlib.pyplot as plt

BZ = np.array(cfg.brillouin_zone)

def plot_gap_mesh(dataset):
    x, y = list(), list()
    gap = firefly.field.load_field_from_file(dataset[0])
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


