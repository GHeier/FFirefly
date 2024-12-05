import pandas as pd
import numpy as np
#import gapFunctions
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from itertools import chain

def gplot(d,eig_num):
    n = len(d)
    for i in range(n):
        theta = i / n * 3.14159*2
        k = [1,theta]
        k1 = [k[0]*np.cos(k[1]), k[0]*np.sin(k[1])]
        g = phys_funcs.g_function(k)
        ratio = (d[i][0]**2 + d[i][1]**2)**(1/2)
        plt.quiver(*k1, g[0], g[1], color=['b'], scale=16)
        plt.quiver(*k1, d[i][0]/ratio, d[i][1]/ratio, color=['r'], scale=16)
    plt.title('Eig #: ' + str(eig_num+1))
    plt.show()

def plot_DOS(): # Plot DOS
    plt.rcParams.update({'font.size': 18})
    mpl.rcParams['lines.linewidth'] = 2
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42
    file = open("DOS.txt")
    df = pd.read_csv(file, delim_whitespace = True, header=None)
    plt.scatter(df.iloc[:,0], df.iloc[:,1], s=3, c='green')
    plt.xlabel(r"$\mu$")
    plt.ylabel(r"$N(\epsilon)$")
    plt.tight_layout() # to fit everything in the prescribed area
    plt.show()

def plot_eigenvalue_divergence():
    file = open("../tests/eigenvalue_divergence.txt")
    df = pd.read_csv(file, delim_whitespace=True, header=None)
    x1, x2, wc = list(), list(), list()
    error = list()
    for i, row in df.iterrows():
        if i%2 == 0:
            wc.append(df.iloc[i,0])
            x1.append(df.iloc[i,1])
        else:
            x2.append(df.iloc[i,1])

    plt.figure(0)
    plt.xlabel(r"$\omega_C$")
    plt.ylabel(r"$\lambda$")
    plt.plot(wc, x1, label=r"$\omega$ independent")
    plt.plot(wc, x2, label=r"$\omega$ dependent")
    plt.show()

def plot_area():
    df = pd.read_csv("info.log", delim_whitespace=True)
    y = df.iloc[:,0]
    x = range(len(y))
    plt.scatter(x,y, s=0.1)
    plt.show()
