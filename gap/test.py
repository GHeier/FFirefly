import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

n = 20
matrix = np.fromfunction(lambda i,j: 0.0+np.cos(i/n*6.28)*np.cos(j/n*6.28), (n,n), dtype=float)
w, v = np.linalg.eig(matrix)
max = 0
max2 = 0
vec = v[:,0].real
vec2 = v[:,0].real
for i in range(len(w)):
    if max < w[i]:
        max2 = max
        max = w[i]
        vec2 = vec
        vec = v[:,i].real
print(max, vec)
sum = 0
for i in range(len(vec)):
    sum += np.cos(i/n*6.28)**2
print(sum)
#print(max2, vec2)
plt.plot(range(n), vec)
#plt.plot(range(n), vec2)
plt.show()
