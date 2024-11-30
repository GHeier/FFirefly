import numpy as np
import matplotlib.pyplot as plt

def test():
    g = np.array([0,1,0])
    for i in range(100):
        for j in range(1,100):
            for k in range(100):
                    Q = get_Q_matrix(i/10,j/10,k/100,g)
                    eigs, vecs = np.linalg.eig(Q)
                    a, b = a__b(i/10,j/10,k/100,g)
                    #print(a/min(eigs))

def test2():
    a_list, b_list = list(), list()
    count, countBigger = 0,0
    for i in range(2):
        theta = np.pi * i / 10
        g = np.array([np.cos(theta),np.sin(theta),0])
        for j in range(1,100):
            T = j / 10
            for k in range(500):
                for l in range(10):
                    alpha = l / 10
                    a, b = a__b( (k/100)**2, T, alpha, g)
                    a_list.append(a)
                    b_list.append(b)
                    s = b - a
                    if abs(a) < abs(s):
                        if alpha != 0:
                            print('hi')
                        #print(theta,T,alpha, k/100)
                        countBigger +=1
                    count +=1
    print(countBigger / count)
    plt.figure()
    plt.plot(a_list)
    plt.plot(b_list)
    plt.show()

