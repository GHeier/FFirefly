import os
import pandas as pd
import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
import numpy as np
import re


to_plot = 'phonon_coulomb' #which files we want to plot, keeping this general in case things change in the future. atm not that important bc we only have one filetype to plot lol



data_directory = '../data/' #move down one level and take from data directory
filenames = [f for f in os.listdir(data_directory) if f.startswith(to_plot)] #look at every file in ../data/ and take only the ones that start with the string specifies by the to_plot variable


#making dictionary that contains the data and the filename for easy call
data_dict = {}
mu_pattern = r"mu=([-\d.]+)" #look for the pattern mu=. then capture everything after such that all following elements are a minus(-), digit (\d), or a period(.). the parenthesis () indicates to capure it; the bracket [] indicates the pattern. The + sign means "one or more." see regex patterns.


mu_vals = []
eig_vals = []

for file in filenames:
    
    # DEEFINING DICT TO ASSOCIATE MU WITH DATA AND WITH FILENAME

    file_path = os.path.join(data_directory, file)
    
    match = re.search(mu_pattern, file) #search each file for the mu pattern
    if match:
        mu_val = float(match.group(1)) #extract the match group. In this case there is only one possible match, but in some cases there are multiple
    else:
        continue # this shouldnt happen, but is kept in case

    data = np.loadtxt(file_path, skiprows=1)
    data_dict[file] = (mu_val, data)
    #the above is kinda unnecessary and should be removed in the future if we still dont use it

    # DEFINING MU VALS AND EIG VALS LIST
    
    eig_vals.append(data[1,3])
    mu_vals.append(mu_val)



#print(len(eig_vals), len(mu_vals))
#print(type(eig_vals), type(mu_vals))
#print(eig_vals, mu_vals)
mu_vals = sorted(mu_vals)
eig_vals = sorted(eig_vals)
plt.figure()
plt.scatter(mu_vals, eig_vals)
plt.xlabel(r"$\mu$")
plt.ylabel(r"$\lambda$")
plt.show()

