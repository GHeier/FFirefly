import numpy as np


def get_path(path, BZ, npts):
    # Convert BZ to numpy array and determine dimension
    BZ = np.array(BZ)
    dim = BZ.shape[0]

    # Define symmetry points (will be trimmed to dimension)
    sympts_3d = {
        'G' : [0.0, 0.0, 0.0],
        'X' : [0.5, 0.0, 0.0],
        'M' : [0.5, 0.5, 0.0],
        'Y' : [0.0, 0.5, 0.0],
        'Z' : [0.0, 0.0, 0.5],
        'R' : [0.0, 0.5, 0.5],
        'A' : [0.5, 0.5, 0.5]
    }

    # Trim to actual dimension
    sympts = {k: v[:dim] for k, v in sympts_3d.items()}

    c1, c2 = "", ""
    portions = []
    for i in range(1, len(path)):
        c1 = path[i - 1].upper()
        c2 = path[i].upper()
        p = np.linspace(sympts[c1], sympts[c2], npts) @ BZ
        portions.append(p)
    return np.vstack(portions)

def format_path(ax, path, npts):
    ax.set_xticks([])
    for i in range(len(path)):
        # Position: start of each segment, except last point is at end of last segment
        if i == len(path) - 1:
            pos = (len(path) - 1) * npts - 1
        else:
            pos = i * npts

        let = path[i].upper()
        if let == 'G':
            let = 'Γ'
        ax.text(pos, -0.02, let, transform=ax.get_xaxis_transform(), ha="center", va="top")
    ax.set_xlim(0, (len(path) - 1) * npts - 1)
