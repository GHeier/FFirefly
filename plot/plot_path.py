import firefly as fly
import firefly.config as cfg

import numpy as np
import matplotlib
matplotlib.use("TkAgg")
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

def plot_section(field, qi, qf, section):
    N = 100
    #BZ = field.domain
    qi = np.array(qi) / 2.0
    qf = np.array(qf) / 2.0
    #BZ = np.array([[2 * np.pi, 0, 0], [0, 2*np.pi, 0], [0, 0, 2*np.pi]])
    BZ = np.array([[1.6388, 0, 0], [0, 1.615, 0], [0, 0, 0.538]])
    nbnd = 8

    x = np.linspace(0, 1, N)
    y = []

    for n in range(1, nbnd+1):
        temp = []
        for t in x:
            q = qi + t * (qf - qi)   
            q_cart = (BZ @ q).tolist()
            temp.append(field(n, q_cart)) 
            #temp.append(fly.epsilon(1, q_cart)) 
        y.append(temp)

    x = np.linspace(section-1, section, N)
    for n in range(nbnd):
        plt.plot(x, y[n], color='#3887f3')
    plt.axvline(x=section, color='gray', linestyle='-', linewidth=1)
    #plt.plot(x, y, color='#3887f3')
    return np.min(y)

def plot_path(files, path):
    fig, ax = plt.subplots()
    ax.set_xticks([])
    ymin = 0

    for file in files:
        field = fly.Field_R(file)
        letter = path[0].lower()
        qi = BZ_point_to_q(letter)
        i = 1
        while i < len(path):
            letter = path[i].lower()
            qf = BZ_point_to_q(letter)
            ymin = min(ymin, plot_section(field, qi, qf, i))
            qi = qf
            i += 1

    plt.axhline(y=0, color='gray', linestyle='--', linewidth=1)
    plt.xlim(0, len(path) - 1)
    for i in range(len(path)):
        ax.text(i, ymin - 1.0, path[i].upper(), ha='center', va='top')
    #fig.patch.set_facecolor('black') 
    #ax.set_facecolor('black')              # Axes background
    plt.show()

