import os
import matplotlib
matplotlib.use('TKAgg')  # Use a backend compatible with your system
import matplotlib.pyplot as plt
import numpy as np
import re

# Define variables
to_plot = 'phonon_coulomb'
data_directory = '../data/'
filenames = [f for f in os.listdir(data_directory) if f.startswith(to_plot)]

# Regex pattern for extracting mu values
mu_pattern = r"mu=([-\d.]+)"

# Initialize lists for mu values and eigenvalues
mu_vals = []
eig_vals = []

# Loop through files
for file in filenames:
    file_path = os.path.join(data_directory, file)
    match = re.search(mu_pattern, file)
    if match:
        mu_val = float(match.group(1))  # Extract mu value
    else:
        continue

    with open(file_path, "r") as f:
        line = f.readline()
        val = line.split(" ")[3]
    #data = np.loadtxt(file_path, skiprows=0)
    #if data.ndim < 2 or data.shape[1] < 4:
    #    raise ValueError(f"Data in {file_path} must have at least two rows and four columns.")

    #eig_vals.append(data[0, 3])  # Extract the eigenvalue
        eig_vals.append(float(val))  # Extract the eigenvalue
        mu_vals.append(mu_val)

# Ensure the relationship between mu_vals and eig_vals is preserved
mu_vals, eig_vals = zip(*sorted(zip(mu_vals, eig_vals)))
print(eig_vals)
# Plot the data
plt.figure()
plt.scatter(mu_vals, eig_vals)
plt.xlabel(r"$\mu$")
plt.ylabel(r"$\lambda$")
plt.show()

