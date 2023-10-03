import pandas as pd
import numpy as np
#import gapFunctions
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from itertools import chain
#import plotly.graph_objects as go

# Compare d-wave direction with g-vector direction
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

# Plot gap function as a function of angle only
def plotGapFunction(eig, gap, num):
    theta = [i*2*np.pi/len(gap) for i in range(0, len(gap))]
    g = [val for val in gap]
    plt.title("Gap Function(Theta) #"+str(num)+", Eig="+str(round(eig,3)))
    plt.xlabel("Angle")
    plt.ylabel("Gap Function")
    plt.plot(theta, g)
    plt.show()

def spherical_to_cartesian(r, theta, phi):
    x = r*np.cos(theta)*np.sin(phi)
    y = r*np.sin(theta)*np.sin(phi)
    z = r*np.cos(phi)
    return x, y, z

def plotGap(eigNum, file):
    df = pd.read_csv(file, delim_whitespace=True)

    if 'phi' in df.columns:
        plt.plot(df.loc[:,'theta'], df.loc['phi'], df.iloc[:,3+eigNum])
    else:
        plt.plot(df.loc[:,'theta'], df.iloc[:,2+eigNum])
        plt.title('$\Delta(\\theta, \phi)$')
        plt.xlabel('$\\theta$')

    plt.ylabel('$\Delta$')
    plt.show()

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
    plt.title("Plotting Scalapino potential")
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

def plot_phi_gap(eigNum, file):
    df = pd.read_csv(file, delim_whitespace=True)

    theta_vals = df.theta.unique()
    for i in range(len(theta_vals)):
        phis, values = list(), list()
        for j in range(len(df.loc[:,'phi'])):
            x = theta_vals[i]
            if x == df.loc[j, 'theta']:
                phis.append(df.loc[j,'phi'])
                values.append(df.iloc[j,3+eigNum])

        #plt.plot(thetas, values, marker='o', ms=1, linestyle='None')
        plt.plot(phis, values)

        plt.title('$\Delta(\\phi), \\theta=$'+str(x))
        plt.xlabel('$\\phi$')

        plt.ylabel('$\Delta$')
        plt.show()

def plot_theta_gap(eigNum, file):
    df = pd.read_csv(file, delim_whitespace=True)

    phi_vals = df.phi.unique()
    for i in range(len(phi_vals)):
        thetas, values = list(), list()
        plt.figure(i)
        for j in range(len(df.loc[:,'theta'])):
            x = phi_vals[i]
            if x == df.loc[j, 'phi']:
                thetas.append(df.loc[j,'theta'])
                values.append(df.iloc[j,3+eigNum])

        #plt.plot(thetas, values, marker='o', ms=1, linestyle='None')
        plt.plot(thetas, values)

        plt.title('$\Delta(\\theta), \\phi=$'+str(x))
        plt.xlabel('$\\theta$')

        plt.ylabel('$\Delta$')
        plt.show()

def plotFS(file):
    df = pd.read_csv(file, delim_whitespace=True)

    fig = plt.figure()
    if 'phi' in df.columns:
        ax = fig.add_subplot(projection='3d')
        for i in range(len(df)):
            r, theta, phi = df.loc[i,'r'], df.loc[i, 'theta'], df.loc[i, 'phi']
            x, y, z = spherical_to_cartesian(r, theta, phi)
            x, y, z = r, theta, phi;
            ax.scatter(x, y, z, s=1, c='green')
    else:
        ax = fig.add_subplot(projection='polar')
        ax.yaxis.get_major_locator().base.set_params(nbins=5)
        ax.spines['polar'].set_visible(False)
        #ax.grid(False)
        plt.polar(df.loc[:,'theta'], df.loc[:,'r'], "o", markersize=1)

    plt.title("Fermi Surface")
    plt.show()

def test(file):
    df = pd.read_csv(file, delim_whitespace=True)

    fig = plt.figure()
    if 'phi' in df.columns:
        ax = fig.add_subplot(projection='3d')
        for i in range(len(df)):
            if abs(df.loc[i, 'phi']) > 0.05:
                continue
            r, theta, phi = df.loc[i,'r'], df.loc[i, 'theta'], df.loc[i, 'phi']
            x, y, z = spherical_to_cartesian(r, theta, phi)
            ax.scatter(x, y, z, s=1, c='green')
    else:
        ax = fig.add_subplot(projection='polar')
        ax.yaxis.get_major_locator().base.set_params(nbins=5)
        ax.spines['polar'].set_visible(False)
        #ax.grid(False)
        plt.polar(df.loc[:,'theta'], df.loc[:,'r'], "o", markersize=1)

    plt.title("Fermi Surface")
    plt.show()

def plot_chi_test(file):
    df = pd.read_csv(file, delim_whitespace=True)
    size = int(len(df)/4)
    q_vals = df.iloc[0:size,0]
    plt.figure(1)
    plt.title("Chi(q) @ T=750K along (1,1,1)")
    for i in range(4):
        col = df.iloc[i*size:(i+1)*size,1]
        plt.plot(q_vals, col)
    plt.legend(["mu=0", "mu=-1", "mu=-2", "mu=-3"])
    plt.ylabel(r'$\chi$')
    plt.xlabel(r'q')
    plt.show()

def plot_chi_single_test(file):
    df = pd.read_csv(file, delim_whitespace=True)
    q_vals = df.iloc[:,0]
    plt.figure(1)
    col = df.iloc[:,1]
    plt.scatter(q_vals, col, s=0.5)
    plt.show()


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

def plotGap_cart(potential, n, mu, dim, U, eigNum):
    file = "../data/"+str(potential)+str(dim)+"D_mu="+str(mu)+"_U="+str(U)+"_wD=1.0_n="+str(n)+".dat"
    df = pd.read_csv(file, delim_whitespace=True)
    fig = plt.figure(0)
    ax = fig.add_subplot(111, projection='3d')
    x = df.loc[:,'x'] 
    y = df.loc[:,'y'] 
    #z = df.loc[:,'z'] 
    c = df.iloc[:,2+eigNum]
    img = ax.scatter(x, y, c)
    plt.title('$\Delta$ @ $\mu=$'+str(mu))
    plt.xlabel('kx')
    plt.ylabel('ky')
    plt.show()

def plot_4D_gap(file):
    df = pd.read_csv(file, delim_whitespace=True)

    x = df.iloc[:,0]
    y = df.iloc[:,1]
    z = df.iloc[:,2]

    for i in range(2):
        c = df.iloc[:,3+i]
        fig = plt.figure(0)
        ax = fig.add_subplot(111, projection='3d')
        img = ax.scatter(x, y, z, c=c, cmap=cm.coolwarm) #mpl.colormaps['plasma'])
        fig.colorbar(img)
        plt.title(r"$\lambda$ #" + str(i) + r"$\Delta$ @ $\mu=-1.2$") 
        plt.xlabel(r'$k_x$')
        plt.ylabel(r'$k_y$')
        ax.set_zlabel(r'$k_z$')
        plt.show()

def plot_coupling(file):
    df = pd.read_csv(file, delim_whitespace = True, header=None)
    mu = df.iloc[:,0]
    renormalization = df.iloc[:,1]
    dx2y2 = df.iloc[:,2]
    dz2_1 = df.iloc[:,3]
    dxy = df.iloc[:,4]
    dxz = df.iloc[:,5]
    dyz = df.iloc[:,6]
    px = df.iloc[:,7]
    py = df.iloc[:,8]
    pz = df.iloc[:,9]
    s = df.iloc[:,10]
    ext_s = df.iloc[:,11]
    eig = df.iloc[:,12]
    #dx = df.iloc[:,10]
    #dy = df.iloc[:,11]
    #dz = df.iloc[:,12]
    plt.plot(mu, dx2y2, label='dx2y2')
    plt.plot(mu, dz2_1, label='dz2-1')
    plt.plot(mu, dxy , label='dxy')
    plt.plot(mu, dxz , label='dxz')
    plt.plot(mu, dyz , label='dyz')
    plt.plot(mu, px , label = 'px')
    plt.plot(mu, py , label='py')
    plt.plot(mu, s , label='s')
    plt.plot(mu, ext_s , label='ext_s')
    plt.plot(mu, eig , label='eig')
    #plt.plot(mu, dx , label = 'dx')
    #plt.plot(mu, dy )
    #plt.plot(mu, dz )
    plt.title("Coupling Constants")
    plt.xlabel('mu')
    plt.ylabel('$\lambda$')
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

def plot_chi_temp(file): 
    df = pd.read_csv(file, delim_whitespace=True)
    x = df.iloc[:,0]
    y = df.iloc[:,1]
    #z = df.iloc[:,2]
    plt.scatter(x, y, s=1.0)
    #plt.scatter(x, z, s=10.0)

def plot_1D_test(file):
    df = pd.read_csv(file, delim_whitespace=True, header=None)
    x = df.iloc[:,0]
    y = df.iloc[:,1]
    plt.scatter(x, y, s=10.0)
    plt.show()

def plot_ep(file):
    df = pd.read_csv(file, delim_whitespace=True, header=None)
    x = df.iloc[:,0]
    plt.plot(x)
    plt.show()

def plot_area():
    df = pd.read_csv("info.log", delim_whitespace=True)
    y = df.iloc[:,0]
    x = range(len(y))
    plt.scatter(x,y, s=0.1)
    plt.show()

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

class MyStruct():
    def __init__(self, field1, field2, field3):
        self.q = field1
        self.chi = field2
        self.category = field3

def get_SAV_chi(file, orbital):
    f = open(file, "r")
    q = list()
    values = list()
    last_category = ""
    for x in f:
        if ';' in x:
            last_category = x
        if "q-vector" in x:
            x = x.split()
            for i in range(3):
                q.append(float(x[i]))

        if x.count(orbital) == 4:
            x = x.split()
            if x[0] == "(":
                values.append(MyStruct(q, float(x[1][:-1]), last_category))
            else:
                values.append(MyStruct(q, float(x[0][1:-1]), last_category))
    return values

def get_categories(file):
    f = open(file, "r")
    categories = list()
    for x in f:
        if ';' in x:
            categories.append(x)
    return categories

def is_float(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

def plot_trace_SAV_chis():
    file_base = "/home/g/PhDResearch/GapFunction/LmtART913_3/CHI/chi"
    chis = list()
    for i in range(1, 351):
        file = file_base
        if i < 100:
            file = file_base + "0"
        if i < 10:
            file += "0"
        f = open(file+str(i), "r")
        sum = 0
        for x in f:
            x = x.split()
            for string in x:
                if is_float(string):
                    sum += float(string)
        chis.append(sum)
    plt.plot(chis)
    plt.show()

def plot_all_SAV_chis():
    file_base = "/home/g/PhDResearch/GapFunction/LmtART913_2/CHI/chi0"
    orbitals = ["3z2-1", "x2-y2"]
    categories = get_categories(file_base+"01")
    colors = [[1, 0, 0], [0, 1, 0], [0, 0, 1], [0.5, 0.5, 0]]
    color_map = lambda x : [abs(np.cos(x)), abs(np.sin(x)), abs(np.cos(x+1)) ]
    for orbital in orbitals:
        fig = plt.figure(0)
        ax = fig.add_subplot(111, projection='3d')
        plt.title("SAV Chi(q) for " + orbital)
        x, y, z, c = list(), list(), list(), list()# np.zeros((65,3)) 
        for i in range(1, 66):
            file = file_base
            if i < 10:
                file = file_base + "0"
            values = get_SAV_chi(file+str(i), orbital)
            for val in values:
                x.append(val.q[0])
                y.append(val.q[1])
                z.append(val.chi)
                index = categories.index(val.category) 
                c.append(color_map(index))
            #for val in values:
            #    x.append(q[1])
            #    y.append(q[0])
            #    z.append(val)
        img = ax.scatter(x,y,z)
        #img = ax.scatter(0,0,0,color=[1,0,0])
        plt.show()

if __name__ == "__main__":
    #plot_all_SAV_chis()
    plot_trace_SAV_chis()
    #plot_V_symmetry_line()
    #plot_V_terms()
    #plot_coupling("coupling.dat")
    #plot_chi_test("chi_plot2.dat")
    #plot_coupling("coupling.dat")
    #plotGap_cart("scalapino", 20, 0.0, 2, 1.0, 0)
    #plotGap_cart("scalapino", 20, -0.1, 2, 1.0, 0)
    #plotGap_cart("scalapino", 20, -0.2, 2, 1.0, 0)
    #plotGap_cart("scalapino", 20, -0.3, 2, 1.0, 0)
    #plotGap_cart("scalapino", 20, -0.4, 2, 1.0, 0)
    #plotGap_cart("scalapino", 40, -0.5, 2, 4.0, 0)
    #plotGap_cart("scalapino", 40, -0.9, 2, 4.0, 0)
    #plotGap_cart("scalapino", 40, -1.0, 2, 4.0, 0)
    #plotGap_cart("scalapino", 40, -1.1, 2, 4.0, 0)
    #plotGap_cart("scalapino", 40, -1.2, 2, 4.0, 0)
    #plotGap_cart("scalapino", 40, -1.3, 2, 4.0, 0)
    #plotGap_cart("scalapino", 40, -1.4, 2, 4.0, 0)
    #plotGap_cart("scalapino", 40, -1.5, 2, 4.0, 0)
    #plotGap_cart("scalapino", 40, -1.6, 2, 4.0, 0)
    #plotGap_cart("scalapino", 40, -1.7, 2, 4.0, 0)
    #plotGap_cart("scalapino", 40, -2.2, 2, 4.0, 0)
    #plotGap_cart("test", 40, -1.0, 2, 4.0, 0)
