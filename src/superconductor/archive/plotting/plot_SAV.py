import pandas as pd
import numpy as np
#import gapFunctions
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from itertools import chain

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

def struct_to_df(struct):
    x, y, z = list(), list(), list()
    chi = list()
    category = list()
    for val in struct:
        x.append(val.q[0])
        y.append(val.q[1])
        chi.append(val.chi)
        category.append(val.category)
    df = pd.DataFrame({'x':x, 'y':y, 'chi':chi, 'category':category})
    return df

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

def plot_all_SAV_chis(file_addon, title):
    #file_base = "/home/g/PhDResearch/LmtART/LmtART_newest/CHI/chi0"
    file_base = "/home/g/PhDResearch/bcs_diagonalization/data/Nickelate/"+file_addon+"/hii"
    orbitals = ["x2-y2"]
    #categories = get_categories(file_base+"01")
    #colors = [[1, 0, 0], [0, 1, 0], [0, 0, 1], [0.5, 0.5, 0]]
    #color_map = lambda x : [abs(np.cos(x)), abs(np.sin(x)), abs(np.cos(x+1)) ]
    for orbital in orbitals:
        x, y, z, c = list(), list(), list(), list()# np.zeros((65,3)) 
        for i in range(1, 198):
            file = file_base
            if i < 10:
                file = file_base + "0"
            if i < 100:
                file += "0"
            values = get_SAV_chi(file+str(i), orbital)
            for val in values:
                x.append(val.q[0])
                y.append(val.q[1])
                z.append(val.chi)
                #index = categories.index(val.category) 
                #c.append(color_map(index))
            #for val in values:
            #    x.append(q[1])
            #    y.append(q[0])
            #    z.append(val)
        fig = plt.figure(0)
        ax = fig.add_subplot(111, projection='3d')
        plt.title(title)
        plt.xlabel(r'$k_x$')
        plt.ylabel(r'$k_y$')
        img = ax.scatter(x,y,z)
        #img = ax.scatter(0,0,0,color=[1,0,0])
        plt.show()

def symmetry_line_condition(val, line):
    if line == "gamma_m":
        return abs(val.q[0] - val.q[1]) < 0.01
    if line == "x_m":
        return abs(val.q[0] - 0.5) < 0.01
    if line == "gamma_x":
        return abs(val.q[1] - 0.00) < 0.01

def symmetry_line_result(val, line):
    if line == "gamma_x":
        return val.q[0] + 1.0
    elif line == "x_m":
        return 1 - val.q[1]
    else:
        return val.q[0]

def plot_symmetry_SAV_chis(file_addon, title):
    orbitals = ["x2-y2"]
    file_base = "/home/g/PhDResearch/bcs_diagonalization/data/HgBa2CuO4-static/"+file_addon+"/hii"
    #file_base = "/home/g/PhDResearch/LmtART/LmtART_newest/CHI/chi0"

    fig = plt.figure(0)
    plt.title(title)
    spots = [0, 0.5, 1.0, 1.5]
    labels = ["G", "X", "M", "G"]
    plt.xticks(spots, labels)

    for orbital in orbitals:
        x, y = list(), list()
        for line in ["gamma_m", "x_m", "gamma_x"]:
            for i in range(1, 198):
                file = file_base

                if i < 10:
                    file = file_base + "0"
                if i < 100:
                    file += "0"

                values = get_SAV_chi(file+str(i), orbital)
                for ind in range(len(values)):
                    if (symmetry_line_condition(values[ind], line)):
                        if line == "gamma_x":
                            x.append(values[ind].q[0] + 1.0)
                        elif line == "x_m":
                            x.append(1 - values[ind].q[1])
                        else:
                            x.append(values[ind].q[0])
                        y.append(values[ind].chi)
                        color_string = 'orange'
                        spin = "odd"
                        if values[ind].category[-12:-1] == "Dn-Dn Dn-Dn": 
                            color_string = "red"
                            spin = "even"
                        #elif values[ind].category[-12:-1] == "Up-Up Up-Up": 
                        #    color_string = "blue"
                        #if values[ind].category[-12:-1] == "Up-Up Dn-Dn": 
                        #    color_string = "black"
                        #elif values[ind].category[-12:-1] == "Up-Dn Dn-Dn":
                        #    color_string = "green"
                        #elif values[ind].category[-12:-1] == "Dn-Up Up-Up":
                        #    color_string = "purple"
                        elif values[ind].category[-12:-1] == "Dn-Up Up-Dn":
                            color_string = "blue"
                            spin = "zero"


                        #color_string = 'black'
                        plt.scatter(x[-1], y[-1], s=3.0, color=color_string)

        red_patch = mpatches.Patch(color='red', label='Even Spin')
        blue_patch = mpatches.Patch(color='blue', label='Zero Spin')
        orange_patch = mpatches.Patch(color='orange', label='Odd Spin')
        plt.legend(handles=[red_patch, blue_patch, orange_patch])
        plt.savefig("/home/g/PhDResearch/Figures/susceptibility/symmetry_plots/HgBa2CuO4/" + title + ".png")

def plot_spin_susceptibility(file_addon, title, u):
    orbitals = ["x2-y2"]
    if (u == "15"):
        u = "1.5"
    file_base = "/home/g/PhDResearch/bcs_diagonalization/data/Nickelate/U="+u+"/"+file_addon+"/hii"

    fig = plt.figure(0)
    plt.title(title)
    spots = [0, 0.5, 1.0, 1.5]
    labels = ["G", "X", "M", "G"]
    plt.xticks(spots, labels)
    bilayers = ["cState=Ni@2::3d        ; cState=Ni@1::3d", "cState=Ni@2::3d        ; cState=Ni@2::3d"]

    for orbital in orbitals:
        for layer in bilayers:
            even, zero, x1, x2 = list(), list(), list(), list()
            for line in ["gamma_m", "x_m", "gamma_x"]:
                for i in range(1, 351):
                    file = file_base

                    if i < 10:
                        file = file_base + "0"
                    if i < 100:
                        file += "0"

                    values = get_SAV_chi(file+str(i), orbital)
                    for ind in range(len(values)):
                        if symmetry_line_condition(values[ind], line):
                            if layer in values[ind].category:
                                if values[ind].category[-12:-1] == "Up-Up Up-Up": 
                                    even.append(values[ind].chi)
                                    x1.append(symmetry_line_result(values[ind], line))

                                if values[ind].category[-12:-1] == "Up-Up Dn-Dn": 
                                    zero.append(values[ind].chi)
                                    x2.append(symmetry_line_result(values[ind], line))

            y = list()
            for i in range(len(even)):
                y.append(even[i] - zero[i])
            color = 'blue'
            if layer == bilayers[0]:
                color = 'red'
            plt.scatter(x1, y, s=4.0, color=color)
        red_patch = mpatches.Patch(color='red', label='Different atoms')
        blue_patch = mpatches.Patch(color='blue', label='Same atom')
        plt.legend(handles=[red_patch, blue_patch])
        plt.savefig("/home/g/PhDResearch/Figures/susceptibility/symmetry_plots/Nickelate/" + title + ".png")

def plot_symmetry_SAV_gap():
    file_base = "/home/g/PhDResearch/bcs_diagonalization/gap/dagotto.bcsgap"
    file = open(file_base, 'r')
    for i in range(3):
        x, y, z, gap = list(), list(), list(), list()
        gap_here = False
        for line in file:
            if line[:25] == " BCS Gap Equation for S=0":
                if gap_here:
                    break
                gap_here = True
            if gap_here:
                line = line.split()
                if (is_float(line[0])):
                    kx, ky, kz, delta = float(line[0]), float(line[1]), float(line[2]), float(line[3])
                    if ((i == 0 and abs(kx) < 0.1 and ky >= 0) or (i == 1 and abs(ky - 0.75) < 0.1 and kx >= 0) or (i == 2 and abs(kx - ky) < 0.1)):
                        x.append(kx)
                        y.append(ky)
                        z.append(kz)
                        gap.append(delta)

        fig = plt.figure(0)
        ax = fig.add_subplot(111, projection='3d')
        img = ax.scatter(x, y, gap, s=8) #mpl.colormaps['plasma'])
        if i == 0:
            plt.title(r'$\Gamma$ to $X$')
        if i == 1:
            plt.title(r'$X$ to $M$')
        if i == 2:
            plt.title(r'$M$ to $\Gamma$')
        plt.xlabel(r'$k_x$')
        plt.ylabel(r'$k_y$')
        ax.set_zlabel(r'$\Delta$')
        plt.show()

def plot_2D_SAV_gap(file_addon, title):
    #file_base = "/home/g/PhDResearch/bcs_diagonalization/gap/dagotto.bcsgap"
    file_base = "/home/g/PhDResearch/bcs_diagonalization/data/CaCuO2-static/"+file_addon+".bcsgap"
    #file_base = "/home/g/PhDResearch/bcs_diagonalization/data/HgBa2CuO4-static/"+file_addon+".bcsgap"
    file = open(file_base, 'r')
    x, y, z, gap = list(), list(), list(), list()
    gap_here = False
    for line in file:
        if line[:25] == " BCS Gap Equation for S=0":
            print(line)
            if gap_here:
                break
            gap_here = True
        if gap_here:
            line = line.split()
            if (is_float(line[0])):
                kx, ky, kz, delta = float(line[0]), float(line[1]), float(line[2]), float(line[3])
                x.append(kx)
                y.append(ky)
                z.append(kz)
                gap.append(delta)

    fig = plt.figure(0)
    ax = fig.add_subplot(111, projection='3d')
    img = ax.scatter(x, y, gap, s=8) #mpl.colormaps['plasma'])
    plt.title(title)
    plt.xlabel(r'$k_x$')
    plt.ylabel(r'$k_y$')
    ax.set_zlabel(r'$\Delta$')
    plt.show()

def plot_SAV_gap(file_extra, title):
    #file_base = "/home/g/PhDResearch/bcs_diagonalization/gap/dagotto.bcsgap"
    file_base = "/home/g/PhDResearch/bcs_diagonalization/data/HgBa2CuO4-static/"+file_addon+".bcsgap"
    file = open(file_base, 'r')
    x, y, z, gap = list(), list(), list(), list()
    gap_here = False
    for line in file:
        if line[:25] == " BCS Gap Equation for S=0":
            if gap_here:
                break
            gap_here = True
        if gap_here:
            line = line.split()
            if (is_float(line[0])):
                x.append(float(line[0]))
                y.append(float(line[1]))
                z.append(float(line[2]))
                gap.append(float(line[3]))

    fig = plt.figure(0)
    ax = fig.add_subplot(111, projection='3d')
    img = ax.scatter(x, y, z, c=gap, cmap=cm.coolwarm, s=8) #mpl.colormaps['plasma'])
    fig.colorbar(img)
    plt.title(title)
    plt.xlabel(r'$k_x$')
    plt.ylabel(r'$k_y$')
    ax.set_zlabel(r'$k_z$')
    plt.show()

def plot_test_waves_SAV_gap():
    file_base = "/home/g/PhDResearch/bcs_diagonalization/gap/dagotto.bcsgap"
    file = open(file_base, 'r')
    x, y, z, d_waves, s_waves = list(), list(), list(), list(), list()
    gap_here = False
    for line in file:
        if line[:25] == " BCS Gap Equation for S=0":
            if gap_here:
                break
            gap_here = True
        if gap_here:
            line = line.split()
            if (is_float(line[0])):
                x.append(float(line[0]))
                y.append(float(line[1]))
                z.append(float(line[2]))
                d_wave = np.cos(x[-1]) - np.cos(y[-1])
                s_wave = 1
                if (d_wave > 0):
                    s_wave = -1
                d_waves.append(d_wave)
                s_waves.append(s_wave)

    fig = plt.figure(0)
    ax = fig.add_subplot(111, projection='3d')
    img = ax.scatter(x, y, z, c=d_waves, cmap=cm.coolwarm, s=8) #mpl.colormaps['plasma'])
    fig.colorbar(img)
    plt.title("d-wave")
    plt.xlabel(r'$k_x$')
    plt.ylabel(r'$k_y$')
    ax.set_zlabel(r'$k_z$')
    plt.show()

    fig = plt.figure(0)
    ax = fig.add_subplot(111, projection='3d')
    img = ax.scatter(x, y, z, c=s_waves, cmap=cm.coolwarm, s=8) #mpl.colormaps['plasma'])
    fig.colorbar(img)
    plt.title("s-wave")
    plt.xlabel(r'$k_x$')
    plt.ylabel(r'$k_y$')
    ax.set_zlabel(r'$k_z$')
    plt.show()
