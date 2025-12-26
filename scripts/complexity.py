import numpy as np
import matplotlib.pyplot as plt

n_states = 10
nw = 100
nk = 200
dim = 3
num_objs = 4

def total(n_states, nw, nk, dim, num_objs, ind_dim=4):
    return (n_states**ind_dim * nw * nk ** dim * num_objs) * 10**(-9)

plt.figure()

x = np.linspace(1,10,100)
y = total(x, nw, nk, dim, 1)
plt.plot(x,y,label="wk in 3D")
y = total(x, nw, nk, 2, 1)
plt.plot(x,y,label="wk in 2D")

V = total(x, 1, nk, dim, 1)
S = total(x, 1, nk, dim, 1, 2)
plt.plot(x,V+S,label="k in 3D")


plt.title("Storage Requirements for V(w,k)_{abcd} = w*k*nstates^4")
plt.xlabel("n_states")
plt.ylabel("GB")
plt.ylim(0,100)
plt.legend()
plt.show()
