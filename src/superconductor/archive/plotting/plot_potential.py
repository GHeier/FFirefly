import pandas as pd
import numpy as np
#import gapFunctions
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from itertools import chain

def plot_potential_q(file):
    df = pd.read_csv(file, delim_whitespace=True)
    q = df.iloc[:,0]
    V1 = df.iloc[:,1]
    plt.figure(6)
    plt.title("V(q) @ T=0.25 (Theory)")
    plt.xlabel("p'")
    plt.ylabel("V")
    plt.plot(q, V1)
    mngr = plt.get_current_fig_manager()
    mngr.window.wm_geometry("+1200+800")

def plot_potential_k_k(file):
    df = pd.read_csv(file, delim_whitespace=True)
    plt.figure(5)
    plt.title("V(q) (Actual)")
    plt.ylabel("V")
    plt.xlabel("|q|")
    q_mag = df.iloc[:,0]
    V1 = df.iloc[:,1]
    plt.scatter(q_mag, V1)
    mngr = plt.get_current_fig_manager()
    mngr.window.wm_geometry("+1200+0")

def plot_potential(file):
    df = pd.read_csv(file, delim_whitespace=True)
    plt.title("Potential Matrix")
    for i in range(len(df)):
        row = df.iloc[i]
        plt.plot(range(len(row)), row)
    plt.show()

def plot_V_terms_3D():
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    U = 4.0
    n = 100
    def Vs1(chi):
        return U**2 * chi / (1 - U * chi) 

    def Vs2(chi):
        return U**3 * chi**2 / (1 - U**2 * chi**2) 

    chi_add = np.linspace(0.13, 0.22, n)
    chi_sub = np.linspace(0.13, 0.22, n)

    chi_plus = list()
    for i in range(n):
        chi_plus.append(n*[chi_add[i]])
    chi_plus = list(chain.from_iterable(chi_plus))

    chi_minus = n*(chi_sub.tolist())

    V = list()
    for i in range(n**2):
        v1 = chi_plus[i]
        v2 = chi_minus[i]
        V.append( Vs1(v1) + Vs2(v2) )

    img = ax.scatter(chi_plus, chi_minus, V)
    plt.title("Plotting FLEX potential")
    plt.xlabel("chi_add")
    plt.ylabel("chi_sub")
    plt.show()

def plot_V_terms():
    U = 4.0
    Vs1 = lambda chi : U**2 * chi / (1 - U * chi) 
    Vs2 = lambda chi : U**3 * chi**2 / (1 - U**2 * chi**2) 
    Vt = lambda chi : -U**2 * chi / (1 - U**2 * chi**2)
    chi_vals = np.linspace(0.10, 0.9, 100)
    plt.figure(3)
    plt.title(r"$V(\chi)$ (RPA)")
    plt.xlabel(r"$\chi$")
    plt.ylabel("V")
    #plt.plot(chi_vals, Vs1(chi_vals)+Vs2(chi_vals))
    plt.plot(chi_vals, Vt(chi_vals))
    plt.show()

def plot_V(k_file, V_file):
    k_vals = pd.read_csv(k_file, delim_whitespace=True)
    V_vals = pd.read_csv(V_file, delim_whitespace=True, header=None)

    r = k_vals.iloc[:,0]
    theta = k_vals.iloc[:,1]
    phi = k_vals.iloc[:,2]

    for i in range(10):
        V_col = V_vals.iloc[:,i]
        plt.scatter(r, V_col)
        plt.show()
    
def plot_V_chi(file):
    df = pd.read_csv(file, delim_whitespace=True)
    chi_sub = df.iloc[:,0]
    potential = df.iloc[:,1]
    fig = plt.figure(4)
    plt.title("V(chi) (Actual)")
    plt.scatter(chi_sub, potential)
    mngr = plt.get_current_fig_manager()
    mngr.window.wm_geometry("+600+800")

def plot_V_symmetry_line():
    df = pd.read_csv("info.log", delim_whitespace=True)
    size = int(len(df)/2)
    q_vals = df.iloc[0:size,0]
    plt.figure(1)
    plt.title("Triplet V(q) @ T=2900K along (1,1,0)")
    for i in range(2):
        col = df.iloc[i*size:(i+1)*size,1]
        plt.plot(q_vals, col)
    plt.legend(["mu=0", "mu=-1"])
    plt.ylabel(r'$V(q)$')
    plt.xlabel(r'q')
    plt.show()

