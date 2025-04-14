import ffirefly
import ffirefly.config as cfg

import numpy as np
import matplotlib.pyplot as plt

def BZ_point_to_q(letter):
    if (letter == 'g'): 
        return [0, 0, 0]
    if (letter == 'x'):
        return [1, 0, 0]
    if (letter == 'm'):
        return [1, 1, 0]
    else:
        print("Letter not recognized in BZ")
    return 0

def plot_section(field, qi, qf):
    x, y = list(), list()
    BZ = field.domain
    for i in range(100):
        q = (qf - qi) / 100 + qi
        x.append(BZ * q)
        y.append(x[-1])
    plt.plot(x, y)

def plot_path(file, path):
    fig, ax = plt.subplots()

    field = CMField(file)
    i = 0
    letter = path[i].lower()
    qi = BZ_point_to_q(letter)
    ax.text(qi, -0.1, letter.upper())

    while i < length(path):
        letter = path[i].lower()
        qf = BZ_point_to_q(letter)
        plot_section(field, qi, qf)
        qi = qf
        ax.text(qi, -0.1, letter.upper())

    plt.show()

