import firefly as fly
import firefly.config as cfg

import numpy as np
import os
import matplotlib

# matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

colors=["#9a05fc", "#f00524", "#f0b802", "#3887f3", "#ed6c09", "#20c714", "black", "gray"]

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

    for i in range(len(files)):
        file = files[i]
        label = os.path.splitext(os.path.basename(file))[0]  # gives 'T=0.01'
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
    ax.legend()
    return fig, ax
