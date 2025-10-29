"""
Frequency plotting module for DMFT results.
Plots Green's functions, self-energies, and other frequency-dependent quantities.
"""

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from typing import Optional
import numpy as np
import sys
from pathlib import Path
import firefly as fly


def plot_frequency_data(file_names, ax: Optional[Axes] = None, positive_only=True,
                       plot_real=True, plot_imag=True, labels=None, **kwargs) -> Axes:
    if ax is None:
        fig, ax = plt.subplots()

    for file in file_names:
        field = fly.Field_C(file)
        x = field.w_points
        x = x[x > 0]
        y = field(x)
        y_real = y.real
        y_imag = y.imag

        label = file.split(".h5")[0]
        if np.max(np.abs(y_real)) > 1e-8:
            ax.plot(x, y_real, marker='s', linestyle='-', label=f'Real {label}', markersize=3, linewidth=1.5)
        if np.max(np.abs(y_imag)) > 1e-8:
            ax.plot(x, y_imag, marker='s', linestyle='-', label=f'Imag {label}', markersize=3, linewidth=1.5)
    ax.legend()
    ax.set_xlabel(r'$\omega$')
    return ax
