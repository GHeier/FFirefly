import pandas as pd
import matplotlib.pyplot as plt

def plot_basic(files, **kwargs):
    datasets = load_data(files)
    return line(datasets, **kwargs)

def plot_scatter(files, **kwargs):
    datasets = load_data(files)
    return scatter(datasets, **kwargs)

def plot_hist(files, **kwargs):
    datasets = load_data(files)
    return hist(datasets, **kwargs)

def plot_bar(files, **kwargs):
    datasets = load_data(files)
    return bar(datasets, **kwargs)

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

    return datasets


def line(datasets, **kwargs):
    fig, ax = plt.subplots()

    for x, y, x_label, y_label, label in datasets:
        ax.plot(x, y, label=label, **kwargs)

    ax.set_xlabel(datasets[0][2])
    ax.set_ylabel(datasets[0][3])
    ax.set_title(f"{datasets[0][4]}")
    ax.grid(True)
    ax.legend()

    return fig, ax


def scatter(datasets, **kwargs):
    fig, ax = plt.subplots()

    for x, y, x_label, y_label, label in datasets:
        ax.scatter(x, y, label=label, **kwargs)

    ax.set_xlabel(datasets[0][2])
    ax.set_ylabel(datasets[0][3])
    ax.set_title(f"{datasets[0][4]}")
    ax.grid(True)
    ax.legend()

    return fig, ax


def hist(datasets, **kwargs):
    fig, ax = plt.subplots()

    for x, y, x_label, y_label, label in datasets:
        ax.hist(x, y, label=label, **kwargs)

    ax.set_xlabel(datasets[0][2])
    ax.set_ylabel(datasets[0][3])
    ax.set_title(f"{datasets[0][4]}")
    ax.grid(True)
    ax.legend()

    return fig, ax


def bar(datasets, **kwargs):
    fig, ax = plt.subplots()

    for x, y, x_label, y_label, label in datasets:
        ax.bar(x, y, label=label, **kwargs)

    ax.set_xlabel(datasets[0][2])
    ax.set_ylabel(datasets[0][3])
    ax.set_title(f"{datasets[0][4]}")
    ax.grid(True)
    ax.legend()

    return fig, ax
