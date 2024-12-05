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


def get_file_name(potential, n, mu, dim, U, wc, FS_only):
    file = ("../data/"
            + potential
            + str(dim)
            + "D_mu=" + str(mu)
            + "_U=" + str(U)
            + "_wc=" + str(wc)
            + "_n=" + str(n)
            + "_" + FS_only + ".dat"
            )
    return file

if __name__ == '__main__':
    potential = "phonon_coulomb"
    dim = 2
    n = 20
    mu = -1.0
    U = 4.0
    wc = 0.1
    FS_only = "FS_only"

    file1 = get_file_name(potential, n, mu, dim, U, wc, FS_only)
    print("Plotting ", file1)
    if (dim == 3):
        plot_gap.plot_4D_gap(file1)
    if (dim == 2):
        plot_gap.plotGap_color(file1, mu, 0)
    else:
        print("Invalid dimension")
