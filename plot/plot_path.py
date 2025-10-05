import firefly as fly
import firefly.config as cfg

import numpy as np
import os
import matplotlib
import h5py

# matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

colors=["#9a05fc", "#f00524", "#f0b802", "#3887f3", "#ed6c09", "#20c714", "black", "gray"]

# State colors for eigenvector weights
state_colors = {
    0: "#f00524",   # Red for first state (1,0,0)
    1: "#3887f3",   # Blue for second state (0,1,0)
    2: "#20c714"    # Green for third state (0,0,1)
}

def load_and_check_matrix_field(filename):
    """Load field and check if it's a matrix field."""
    with h5py.File(filename, 'r') as f:
        n_indices = f['n_indices'][()]
        dim_indices = f['dim_indices'][()]
        is_matrix = (n_indices == 2)

        if is_matrix:
            # Load matrix data
            mesh = f['mesh'][:]
            domain = f['domain'][:]
            real_data = f['values/real'][:]
            is_complex = f['is_complex'][()]

            if is_complex:
                imag_data = f['values/imag'][:]
                data = real_data + 1j * imag_data
            else:
                data = real_data

            # Reshape to [num_matrices, mat_dim, mat_dim]
            num_matrices = len(data) // (dim_indices * dim_indices)
            data = data.reshape(num_matrices, dim_indices, dim_indices)

            return {
                'is_matrix': True,
                'data': data,
                'mesh': mesh,
                'domain': domain,
                'mat_dim': dim_indices,
                'num_matrices': num_matrices
            }
        else:
            return {'is_matrix': False}

def diagonalize_along_path(matrix_data, qi, qf, N, BZ):
    """Diagonalize matrices along a path and return eigenvalues and eigenvectors."""
    mesh = matrix_data['mesh']
    mat_dim = matrix_data['mat_dim']
    nx, ny, nz = mesh

    eigenvalues = []
    eigenvectors = []

    for t in np.linspace(0, 1, N):
        q = qi + t * (qf - qi)
        q_cart = BZ @ q

        # Map to mesh indices
        ix = int(q[0] * nx) % nx
        iy = int(q[1] * ny) % ny
        iz = int(q[2] * nz) % nz

        # Get matrix index
        mat_idx = ix * ny * nz + iy * nz + iz

        # Get matrix
        H = matrix_data['data'][mat_idx]

        # Diagonalize
        evals, evecs = np.linalg.eigh(H)
        eigenvalues.append(evals)
        eigenvectors.append(evecs)

    return np.array(eigenvalues), np.array(eigenvectors)

def get_band_color(eigenvector):
    """Get color based on eigenvector weight."""
    # eigenvector is complex, get absolute values
    weights = np.abs(eigenvector)**2

    # Mix colors based on weights
    color = np.zeros(3)
    hex_to_rgb = lambda h: np.array([int(h[1:3], 16), int(h[3:5], 16), int(h[5:7], 16)]) / 255.0

    for i, w in enumerate(weights):
        if i < 3:  # Only first 3 states
            rgb = hex_to_rgb(state_colors[i])
            color += w * rgb

    return tuple(color)

def BZ_point_to_q(letter):
    if letter == "g":
        return [0, 0, 0]
    if letter == "x":
        return [1, 0, 0]
    if letter == "m":
        return [1, 1, 0]
    if letter == "r":
        return [1, 1, 1]
    else:
        print("Letter not recognized in BZ")
    return 0


def plot_section(field, qi, qf, section, letter, multicolor, color, label=None):
    N = 100
    # BZ = field.domain
    qi = np.array(qi) / 2.0
    qf = np.array(qf) / 2.0
    BZ = np.array([[2 * np.pi, 0, 0], [0, 2*np.pi, 0], [0, 0, 2*np.pi]])
    #BZ = np.array([[1.6388, 0, 0], [0, 1.615, 0], [0, 0, 0.538]])
    nbnd = 1

    x = np.linspace(0, 1, N)
    y = []
    y_width = []

    with_n = False
    for n in range(1, nbnd + 1):
        temp = []
        widths = []
        for t in x:
            q = qi + t * (qf - qi)
            q_cart = (BZ @ q).tolist()
            if with_n:
                val = field(n, q_cart)
            else:
                val = field(q_cart)
            width = val.imag
            val = val.real
            temp.append(val)
            widths.append(width)
        y.append(temp)
        y_width.append(widths)
        if with_n:
            break

    x = np.linspace(section - 1, section, N)
    for n in range(nbnd):
        py = np.array(y[n])
        pw = np.array(y_width[n])
        print(np.max(pw))
        if multicolor:
            plt.plot(x, py, color=color, label=label)
            plt.fill_between(x, py - pw, py + pw,
                            color=color, alpha=0.3)
        else:
            plt.plot(x, y[n], color="#9a05fc")
            plt.fill_between(x, py - pw, py + pw,
                            color="#9a05fc", alpha=0.3)
    plt.axvline(x=section, color="gray", linestyle="-", linewidth=1)
    # plt.plot(x, y, color='#3887f3')
    return np.min(y)


def plot_path(files, path, hline=False, multicolor=False):
    fig, ax = plt.subplots()
    ax.set_xticks([])

    BZ = np.array([[2 * np.pi, 0, 0], [0, 2*np.pi, 0], [0, 0, 2*np.pi]])

    for i in range(len(files)):
        file = files[i]
        label = os.path.splitext(os.path.basename(file))[0]

        # Check if file is matrix field
        matrix_data = load_and_check_matrix_field(file)

        if matrix_data['is_matrix']:
            # Plot eigenvalues with eigenvector coloring
            letter = path[0].lower()
            qi = np.array(BZ_point_to_q(letter)) / 2.0
            section = 1

            for j in range(1, len(path)):
                letter = path[j].lower()
                qf = np.array(BZ_point_to_q(letter)) / 2.0

                N = 100
                evals, evecs = diagonalize_along_path(matrix_data, qi, qf, N, BZ)
                x = np.linspace(section - 1, section, N)

                # Plot each band with color based on eigenvector
                for band in range(matrix_data['mat_dim']):
                    for k in range(len(x)):
                        if k < len(x) - 1:
                            # Get color from eigenvector
                            evec = evecs[k, :, band]
                            color = get_band_color(evec)

                            plt.plot(x[k:k+2], evals[k:k+2, band], color=color, linewidth=2)

                plt.axvline(x=section, color="gray", linestyle="-", linewidth=1)
                qi = qf
                section += 1

        else:
            # Original Field_C plotting
            field = fly.Field_C(file)
            letter = path[0].lower()
            qi = BZ_point_to_q(letter)
            j = 1
            while j < len(path):
                letter = path[j].lower()
                qf = BZ_point_to_q(letter)
                plot_section(field, qi, qf, j, letter, multicolor, colors[i], label=label)
                label = None
                qi = qf
                j += 1

    if hline:
        plt.axhline(y=0, color="gray", linestyle="--", linewidth=1)
    plt.xlim(0, len(path) - 1)
    yloc = -0.02
    for i in range(len(path)):
        let = path[i].upper()
        if let == 'G':
            let = 'Γ'
        ax.text(i, yloc, let, transform=ax.get_xaxis_transform(), ha="center", va="top")
    # fig.patch.set_facecolor('black')
    # ax.set_facecolor('black')              # Axes background
    if not matrix_data.get('is_matrix', False):
        ax.legend()
    return fig, ax
