import firefly
import matplotlib.pyplot as plt
import sys

file = sys.argv[1]

fig, ax = firefly.sketch([file])
plt.show()

