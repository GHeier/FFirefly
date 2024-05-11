import sys
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter



def default_plot(labels=['x','y1','y2','y3','y4'], axes=['x_axis', 'y_axis']):
    # Things to edit for plot
    plt.figure(figsize=(6.4,4.8))
    x_axis, y_axis = axes 
    line_labels = labels 

    file = open(sys.argv[1])
    df = pd.read_csv(file, delim_whitespace=True, header=None)

    plt.rcParams.update({'font.size': 18})
    mpl.rcParams['lines.linewidth'] = 2
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42
    num_axis_points = 4

    ax = plt.axes()

    plt.xlabel(x_axis)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.set_xticks(np.linspace(df.iloc[0,0],df.iloc[-1,0],num_axis_points))

    plt.ylabel(y_axis)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.set_yticks(np.linspace(min(df.min()[1:]),max(df.max()[1:]),num_axis_points))

    for i in range(len(df.columns)-1):
        ax.plot(df.iloc[:,0], df.iloc[:,i+1])
        ax.annotate(line_labels[i], 
                    xy = (df.iloc[-1,0], df.iloc[-1,i+1]), 
                    xytext = (1.02 * df.iloc[-1,0], df.iloc[-1,i+1]) 
                    )

    #plt.legend(ncol=2, loc=1) # 9 means top center
    plt.tight_layout() # to fit everything in the prescribed area
    plt.show()
    save_loc = '/home/g/Research/Papers/Griffin/'
    #plt.savefig(save_loc+'eigenvalue_divergence.png', dpi=300)


if __name__ == "__main__":
    default_plot()
