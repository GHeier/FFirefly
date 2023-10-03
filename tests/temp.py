import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(1, '/home/g/Research/GapFunction/gap')
import create_eq

def eq(x):
    T = 1
    return np.tanh(x/T)

k = np.array([1,2])
print(create_eq.spherical_to_cartesian(k))
