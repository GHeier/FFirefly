import pandas as pd
import numpy as np
#import gapFunctions
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from itertools import chain

def plot_coupling(file):
    df = pd.read_csv(file, delim_whitespace = True, header=None)
    mu = df.iloc[:,0]
    renormalization = df.iloc[:,1]
    dx2y2 = df.iloc[:,2]
    #dz2_1 = df.iloc[:,3]
    dxy = df.iloc[:,3]
    #dxz = df.iloc[:,5]
    #dyz = df.iloc[:,6]
    px = df.iloc[:,4]
    #py = df.iloc[:,8]
    #pz = df.iloc[:,9]
    #s = df.iloc[:,10]
    ext_s = df.iloc[:,5]
    eig = df.iloc[:,6]
    #dx = df.iloc[:,10]
    #dy = df.iloc[:,11]
    #dz = df.iloc[:,12]
    plt.plot(mu, dx2y2, label='dx2y2')
    #plt.plot(mu, dz2_1, label='dz2-1')
    plt.plot(mu, dxy , label='dxy')
    #plt.plot(mu, dxz , label='dxz')
    #plt.plot(mu, dyz , label='dyz')
    plt.plot(mu, px , label = 'px')
    #plt.plot(mu, py , label='py')
    #plt.plot(mu, s , label='s')
    plt.plot(mu, ext_s , label='ext_s')
    plt.plot(mu, eig , label='eig')
    #plt.plot(mu, dx , label = 'dx')
    #plt.plot(mu, dy )
    #plt.plot(mu, dz )
    plt.title("Coupling Constants")
    plt.xlabel(r'$\mu$')
    plt.ylabel(r'$\lambda$')
    plt.xlim(-4.8,-0.4)
    plt.ylim(-0.05,0.10)
    plt.legend()
    plt.show()

def plot_coupling_manual():
    mu, d1, d2, s_ext = list(), list(), list(), list()
    mu = -1.2, -1.4, -1.6, -1.8, -2.0, -2.2, -2.4, -2.6, -2.8, -3.0, -3.2
    d1 = 0.142, 0.106, 0.080, 0.062, 0.044, 0.010, -0.006, -0.020, -0.030, -0.038, -0.041
    d2 = 0.142, 0.106, 0.080, 0.062, 0.044, 0.010, -0.006, -0.020, -0.030, -0.038, -0.041
    s_ext = -1.098, -1.108, -1.116, -1.125, -1.133, -1.138, -1.138, -1.134, -1.127, -1.117, -1.103
    plt.figure()
    plt.plot(mu, d1)
    plt.plot(mu, d2)
    #plt.plot(mu, s_ext)
    plt.show()

