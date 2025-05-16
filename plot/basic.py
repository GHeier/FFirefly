import pandas as pd
import matplotlib.pyplot as plt

def plot_basic(files):
    datasets = load_data(files)
    line(datasets)

def plot_scatter(files):
    datasets = load_data(files)
    scatter(datasets)


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
        x_label = columns[0]
        y_label = columns[1]
        x = df.iloc[:, 0]
        y = df.iloc[:, 1]
        title = clean_filename(file)

        datasets.append((x, y, x_label, y_label, title))

    return datasets


def line(datasets):
    """
    Plot each dataset as a line plot.
    """
    plt.figure()
    for x, y, x_label, y_label, label in datasets:
        plt.plot(x, y, label=label)
    plt.xlabel(datasets[0][2])
    plt.ylabel(datasets[0][3])
    plt.title(f"{datasets[0][4]}")
    plt.grid(True)
    plt.legend()
    plt.show()


def scatter(datasets):
    """
    Plot each dataset as a scatter plot.
    """
    plt.figure()
    for x, y, x_label, y_label, label in datasets:
        plt.scatter(x, y, label=label)
        #plt.scatter(x, y, label=label, s = 20, marker='.')
    plt.xlabel(datasets[0][2])
    plt.ylabel(datasets[0][3])
    plt.title(f"{datasets[0][4]}")
    plt.grid(True)
    plt.legend()
    plt.show()
