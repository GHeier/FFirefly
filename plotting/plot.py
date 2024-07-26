import pandas as pd
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from itertools import chain

import plot_SAV
import plot_gap
import plot_potential
import plot_info
import plot_chi


def get_file_name(potential, n, mu, dim, U, wc):
    file = ("../data/"
            + potential
            + str(dim)
            + "D_mu=" + str(mu)
            + "_U=" + str(U)
            + "_wc=" + str(wc)
            + "_n=" + str(n)
            + ".dat"
            )
    return file

if __name__ == '__main__':
    potential = "scalapino"
    dim = 3
    n = 19
    mu = -1.0
    U = 4.0
    wc = 0.5

    file1 = get_file_name(potential, n, mu, dim, U, wc)
    print(file1)
    plot_gap.plot_4D_gap(file1)
