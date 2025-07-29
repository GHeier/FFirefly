import matplotlib.pyplot as plt
import sys
from pathlib import Path
import re

def find_files_level1(substring, root="."):
    root_path = Path(root)

    # Files in root directory
    root_matches = sorted(root_path.glob(f"*{substring}*"))

    # Files in immediate subdirectories
    subdir_matches = []
    for subdir in sorted([p for p in root_path.iterdir() if p.is_dir()]):
        subdir_matches.extend(sorted(subdir.glob(f"*{substring}*")))

    return root_matches + subdir_matches


def load_and_plot(file_names, **kwargs):
    mu = kwargs.get("mu", None)
    for file_name in file_names:
        with open(file_name, 'r') as f:
            lines = f.readlines()
        title = file_name.split(".")[0]
        title = title.split("/")[-1]

        fermi_energy = mu

        if fermi_energy is None:
            matches = find_files_level1("scf.out")
            if len(matches) == 0:
                fermi_energy = 0.0
            if len(matches) > 0:
                fermifile = matches[0]
                with open(fermifile, 'r') as f:
                    data = f.read()
                    fermi_energy = re.search(r"Fermi energy is\s+([0-9.]+)", data)
                    if fermi_energy:
                        fermi_energy = float(fermi_energy.group(1))
                        print(f"Fermi Energy found at {fermi_energy}")
                    else:
                        print(f"Fermi energy not found in {fermifile}")
                        exit()

        emin = fermi_energy - 4.1
        emax = fermi_energy + 4.1

        x, y = [], []
        
        # Set figure size here (taller than wide)
        fig, ax = plt.subplots(figsize=(4, 8))  # 4 inches wide, 8 inches tall
        
        min_x = 0
        max_x = 0
        for line in lines:
            stripped = line.strip()
            if stripped:
                # Split the line into two parts, x and y
                parts = stripped.split()
                if float(parts[1]) < emin or float(parts[1]) > emax:
                    continue
                x.append(float(parts[0]))
                y.append(float(parts[1]) - fermi_energy)
                if float(parts[0]) < min_x:
                    min_x = float(parts[0])
                if float(parts[0]) > max_x:
                    max_x = float(parts[0])
            else:
                # Empty line, meaning a new segment starts
                ax.plot(x, y, color='red', linewidth=0.8)
                x, y = [], []  # Reset lists for the next segment

        # Plot the last segment if present
        if x and y:
            ax.plot(x, y, color='red', linewidth=0.8)

        # Custom x-axis labels at specific positions with LaTeX for Gamma (Γ)
        xticks = [min_x, min_x + (max_x - min_x) / 3, min_x + 2 * (max_x - min_x) / 3, max_x]
        ax.set_xticks(xticks)  # Four x-tick positions
        ax.set_xticklabels([r'$\Gamma$', 'X', 'M', r'$\Gamma$'])  # LaTeX labels for Gamma, X, M, Gamma

        # Remove x-axis tick marks
        ax.tick_params(axis='x', which='both', length=0)

        # Set plot labels and limits
        ax.set_ylabel('E(k)')
        ax.set_xlim(min_x, max_x)
        ax.set_ylim(-4, 4)
        ax.set_title("YBa2Cu3O7" + " Band Structure")
        ax.grid(True)

    # Display the plot
    plt.show()

# Usage
#file_name = sys.argv[1]
#load_and_plot(file_name, 20)
