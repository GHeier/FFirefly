import pandas as pd
import matplotlib.pyplot as plt
import csv


def plot_basic(files):
    plt.figure()
    for file in files:
        # Step 1: Detect header (look for non-numeric first row)
        with open(file, "r") as f:
            first_line = f.readline().strip()
            try:
                # Try parsing first line as floats — if fails, assume header
                float(first_line.split(None)[0])
                header = None
            except ValueError:
                header = 0

        # Step 2: Read file with pandas
        df = pd.read_csv(file, sep=r"\s+", engine="python", header=header)

        # Step 3: Plot the first two columns
        x = df.iloc[:, 0]
        y = df.iloc[:, 1]

        plt.plot(x, y)
    plt.grid(True)
    plt.show()
