import matplotlib as mpl, scipy as scp, numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from itertools import chain

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def plot_fs(file_path):
    data = np.loadtxt(file_path)

    # data sizing (x, y, z)
    if data.shape[1] != 3:
        raise ValueError("The file must contain 3 columns (x, y, z).")

    x, y, z = data[:, 0], data[:, 1], data[:, 2]

    # 3D scatter plt
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')
    scatter = ax.scatter(x, y, z, c=z, cmap='viridis', s=50)

    # lble
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title("Fermi Surface")

    plt.show()

# Example usage
plot_fs("FS.dat")

