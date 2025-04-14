import os
import sys
import pandas as pd

file = sys.argv[1]
if not os.path.exists(file):
    print("File does not exist")
    sys.exit(1)

num_cols = 2

for i in range(2, len(sys.argv)):
    if sys.argv[i] == "-nl":
        num_cols = int(sys.argv[i + 1])
        break

data = pd.read_csv(file, sep=' ', header=None)

if "bands.dat.gnu" in file:
    pass
elif ".dos" in file:
    pass
else:
    default_plot(num_cols, data)
