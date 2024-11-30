import pandas as pd
import numpy as np
#import gapFunctions
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from itertools import chain

def plot_chi_test(file):
    df = pd.read_csv(file, delim_whitespace=True)
    size = int(len(df)/4)
    q_vals = df.iloc[0:size,0]
    plt.figure(1)
    plt.title(r"$\chi$(q) @ T=300K along (1,1,1)")
    for i in range(4):
        col = df.iloc[i*size:(i+1)*size,1]
        plt.plot(q_vals, col)
    plt.legend(["mu=0", "mu=1", "mu=2", "mu=3"])
    plt.ylabel(r'$\chi$')
    plt.xlabel(r'q')

def plot_chi_single_test(file):
    df = pd.read_csv(file, delim_whitespace=True)
    q_vals = df.iloc[:,0]
    plt.figure(1)
    col = df.iloc[:,1]
    plt.plot(q_vals, col)

def plot_chi_values(file):
    df = pd.read_csv(file, delim_whitespace=True)
    x = df.iloc[:,0]
    y = df.iloc[:,1]
    chi = df.iloc[:,2]
    plt.figure(2)
    plt.title("Chi Scatter plot as function of q (Actual)")
    plt.xlabel("q")
    plt.ylabel("chi")
    plt.scatter(x, chi, s=1)
    plt.show()

def plot_chi_temp(file): 
    df = pd.read_csv(file, delim_whitespace=True)
    x = df.iloc[:,0]
    y = df.iloc[:,1]
    #z = df.iloc[:,2]
    plt.scatter(x, y, s=1.0)
    #plt.scatter(x, z, s=10.0)

