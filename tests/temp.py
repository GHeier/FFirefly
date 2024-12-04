import sys
import numpy as np
import matplotlib.pyplot as plt

def eq1(x):
    return np.cos(x + 0.1) - np.cos(x)

def eq2(x):
    return -2 * np.sin(0.05) * np.sin(x + 0.05)

def eq3(x):
    return 0.0999583*np.cos(x + 1.6208)

def eq4(x):
    return -np.sin(x) * 0.1


x = np.linspace(-3, 3, 1000)
y1 = eq1(x)
y2 = eq2(x)
y3 = eq3(x)
y4 = eq4(x)

plt.figure(0)
plt.plot(x, y1, label='eq1')
plt.plot(x, y2, label='eq2')
plt.plot(x, y3, label='eq3')
plt.plot(x, y4, label='eq4')
plt.legend()
plt.show()
