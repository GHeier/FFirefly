import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from typing import Optional
import firefly as fly

def plot_basic(files, ax: Optional[Axes] = None, **kwargs) -> Axes:
    datasets = load_data(files)
    return line(datasets, ax=ax, **kwargs)

def plot_scatter(files, ax: Optional[Axes] = None, **kwargs) -> Axes:
    datasets = load_data(files)
    return scatter(datasets, ax=ax, **kwargs)

def plot_hist(files, ax: Optional[Axes] = None, **kwargs) -> Axes:
    datasets = load_data(files)
    return hist(datasets, ax=ax, **kwargs)

def plot_bar(files, ax: Optional[Axes] = None, **kwargs) -> Axes:
    datasets = load_data(files)
    return bar(datasets, ax=ax, **kwargs)

def clean_filename(file):
    for ext in [".dat", ".csv", ".txt"]:
        if file.endswith(ext):
            file = file[:-len(ext)]
    return file.replace("_", " ")


def load_data(files):
    """
    Load multiple files, detect headers, and return structured data for plotting.
    Returns: List of tuples (x, y, x_label, y_label, file_name)
    """
    datasets = []

    for file in files:
        # Step 1: Detect header
        if file[-2:] != "h5" and file[-4:] != "hdf5":
            with open(file, "r") as f:
                first_line = f.readline().strip()
                try:
                    float(first_line.split(None)[0])
                    header = None
                except ValueError:
                    header = 0

            # Step 2: Read file
            df = pd.read_csv(file, sep=r"\s+", engine="python", header=header)
            columns = df.columns.tolist()
            if '#' in columns:
                columns.remove('#')
            x_label = columns[0]
            y_label = columns[1]
            x = df.iloc[:, 0]
            y = df.iloc[:, 1]
            title = clean_filename(file)

            datasets.append((x, y, x_label, y_label, title))
        else:
            field = fly.Field_R(file)
            x = field.w_points
            y = field(x)
            datasets.append((x, y, "w", "values", os.path.splitext(os.path.basename(file))[0]))

    return datasets


def line(datasets, ax: Optional[Axes] = None, **kwargs) -> Axes:
    if ax is None:
        fig, ax = plt.subplots()

    for x, y, x_label, y_label, label in datasets:
        ax.plot(x, y, label=label, **kwargs)

    ax.set_xlabel(datasets[0][2])
    ax.set_ylabel(datasets[0][3])
    ax.set_title(f"{datasets[0][4]}")
    ax.grid(True, alpha=0.2)
    ax.legend()

    return ax


def scatter(datasets, ax: Optional[Axes] = None, **kwargs) -> Axes:
    if ax is None:
        fig, ax = plt.subplots()

    for x, y, x_label, y_label, label in datasets:
        ax.scatter(x, y, label=label, **kwargs)

    ax.set_xlabel(datasets[0][2])
    ax.set_ylabel(datasets[0][3])
    ax.set_title(f"{datasets[0][4]}")
    ax.grid(True, alpha=0.2)
    ax.legend()

    return ax


def hist(datasets, ax: Optional[Axes] = None, **kwargs) -> Axes:
    if ax is None:
        fig, ax = plt.subplots()

    for x, y, x_label, y_label, label in datasets:
        ax.hist(x, y, label=label, **kwargs)

    ax.set_xlabel(datasets[0][2])
    ax.set_ylabel(datasets[0][3])
    ax.set_title(f"{datasets[0][4]}")
    ax.grid(True, alpha=0.2)
    ax.legend()

    return ax


def bar(datasets, ax: Optional[Axes] = None, **kwargs) -> Axes:
    if ax is None:
        fig, ax = plt.subplots()

    for x, y, x_label, y_label, label in datasets:
        ax.bar(x, y, label=label, **kwargs)

    ax.set_xlabel(datasets[0][2])
    ax.set_ylabel(datasets[0][3])
    ax.set_title(f"{datasets[0][4]}")
    ax.grid(True, alpha=0.2)
    ax.legend()

    return ax
