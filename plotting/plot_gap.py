import pandas as pd
import numpy as np
#import gapFunctions
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from itertools import chain

# Plot gap function as a function of angle only
def plotGapFunction(eig, gap, num):
    theta = [i*2*np.pi/len(gap) for i in range(0, len(gap))]
    g = [val for val in gap]
    plt.title("Gap Function(Theta) #"+str(num)+", Eig="+str(round(eig,3)))
    plt.xlabel("Angle")
    plt.ylabel("Gap Function")
    plt.plot(theta, g)
    plt.show()

def plotGap(eigNum, file):
    df = pd.read_csv(file, delim_whitespace=True)

    if 'phi' in df.columns:
        plt.plot(df.loc[:,'theta'], df.loc['phi'], df.iloc[:,3+eigNum])
    else:
        plt.plot(df.loc[:,'theta'], df.iloc[:,2+eigNum])
        plt.title(r'$\Delta(\theta, \phi)$')
        plt.xlabel(r'$\theta$')

    plt.ylabel(r'$\Delta$')
    plt.show()

def plotGap_cart(potential, n, mu, dim, U, w_D, eigNum):
    df = pd.read_csv(file, delimiter=' ')
    fig = plt.figure(0)
    ax = fig.add_subplot(111, projection='3d')
    x = df.loc[:,'x'] 
    y = df.loc[:,'y'] 
    #z = df.loc[:,'z'] 
    c = df.iloc[:,2+eigNum]
    img = ax.scatter(x, y, c)
    plt.title(r'$\Delta$ @ $\mu=$'+str(mu))
    plt.xlabel('kx')
    plt.ylabel('ky')
    plt.show()

def plot_4D_gap(file):
    df = pd.read_csv(file, delimiter=' ')

    x = df.iloc[:,0]
    y = df.iloc[:,1]
    z = df.iloc[:,2]

    c = df.iloc[:,3]
    fig = plt.figure(0)
    ax = fig.add_subplot(111, projection='3d')
#        for i in range(len(c)):
#            if abs(c[i]) > 1.0:
#                c[i] = 0.0
    img = ax.scatter(x, y, z, c=c, cmap=cm.coolwarm) #mpl.colormaps['plasma'])
    fig.colorbar(img)
    plt.title(r"$\lambda$ #" + str(0) + r"$\Delta$ @ $\mu=-1.2$") 
    plt.xlabel(r'$k_x$')
    plt.ylabel(r'$k_y$')
    ax.set_zlabel(r'$k_z$')
    plt.show()

