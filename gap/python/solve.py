import numpy as np
import create_eq
import matplotlib.pyplot as plt
from scipy import optimize
from scipy import integrate
from scipy.linalg import eigvals
from scipy.linalg import eig
import pandas as pd
from alive_progress import alive_bar


# angle, radius
angular_points, radial_points, T0, mu = 8, 4, 1, 0
#N = 1


# Sorts eigenvectors based on eigenvalue reordering
# Sorting removes eig-vec pairs with eig=0 and orders by eig size
def sort_eigs_and_vecs(eigs, vecs):
    vecs = np.asarray([vecs[:,i] for i in range(len(eigs)) if np.round(eigs[i], 10) > 0])
    eigs = np.asarray([x for x in eigs if np.round(x,10) > 0])

    indices = sorted(range(len(eigs)), key=lambda k: eigs[k])
    eigs = np.sort(eigs)

    def newVecs(indices, vecs):
        result = np.zeros(shape=vecs.shape)
        for i in range(len(indices)):
            j = indices[i]
            result[len(result)-j-1] = vecs[i].real
        return result.transpose()

    vecs = newVecs(indices, vecs)

    return eigs[::-1].real, vecs 

# Returns:
#   temperature
#   f value
#   eigenvalue list
#   eigenvector list(not transposed)
def solve_simple(V):
    from cfg import numPts, dim
    w, v = np.linalg.eig(V) 
    #print('Eig:\n', w)
    #print('Vec:\n', v)
    w, v = sort_eigs_and_vecs(w, v)

    # Scaling factor
    w *= (1 / numPts) ** dim

    highEig = w[0].real
    f = highEig**(-1)
    Tc = solveTc_energy(f)

    return Tc, w, v

def solve_complex(numPts, dim, f_type, alpha):
    Tc = solveTcComplex(angular_points,radial_points,f_type,alpha)
    P = create_eq.getWorkMatrix(angular_points,radial_points,Tc,f_type,alpha,mu)
    #print (pd.DataFrame(P))
    w, v = np.linalg.eig(P)
    #print('Eig:\n', w)
    #print('Vec:\n', v)
    w, v = sort_eigs_and_vecs(w, v)
    if len(w) > 0 and round(np.max(w),6) != 1:
        exit("Critical Temperature does not exist")

    return [Tc, w, v]

def solve(V,k):
    from cfg import f_type, numPts
    def func(T):
        P = create_eq.get_matrix(V, k, T)
        w = eigvals(P)
        w *= 1/numPts**3
        highEig = np.amax(w).real
        return abs(highEig - 1)
    Tc = optimize.least_squares(func, 50, bounds = (0.01, 1000), verbose=2).x[0] 
    Tc = 1
    P = create_eq.get_matrix(V, k, Tc)
    w, v = eig(P)
    w *= 1/numPts**3
    w, v = sort_eigs_and_vecs(w, v)
    print('Matrix:\n', P)
    print('Eigenvalues:\n', w)
    print('Temperature: ', Tc)
    if len(w) > 0 and round(w[0],6) != 1:
        print("Eig = {}".format(w[0]))
        exit("Critical Temperature does not exist")
    if f_type == 'triplet':
        v.reshape(int(len(v)/dim),dim)
    return Tc, w, v


# Uses optimize package to find the zero of f integral - calculated f value as a function of Tc
def solveTc(invEig):
    from cfg import wc
    def func(Tc):
        def f(x):
            T = Tc[0]
            return np.tanh(x**2/(2*T))/(2*x**2)
        return abs(2*integrate.quadrature(f,0,wc)[0] - invEig)
    return optimize.least_squares(func, 1.0, bounds = (0.01, 1000), xtol = 10**(-14), verbose=0).x[0] 

def solveTc_energy(invEig):
    from cfg import wc
    def func(Tc):
        def f(x):
            T = Tc[0]
            return np.tanh(x/(2*T))/(2*x)
        return abs(2*integrate.quadrature(f,0,wc)[0] - invEig)
    return optimize.least_squares(func, 1.0, bounds = (0.01, 1000), xtol = 10**(-14), verbose=0).x[0] 

def solveTcComplex(angular_points,radial_points,f_type,alpha):
    def func(Tc):
        P = create_eq.getWorkMatrix(angular_points,radial_points,Tc,f_type,alpha,mu)
        #print (pd.DataFrame(P))
        w, v = np.linalg.eig(P)
        highEig = np.amax(w).real
        #print('Eig:  ', highEig)
        return abs(highEig - 1)
    return optimize.least_squares(func, 50, bounds = (0.01, 1000), verbose=0).x[0] 

def save(name, Tc, eigs, gaps, k):
    precision = 3
    table = np.concatenate((k, gaps), axis=1)
    np.savetxt(name, table, fmt='%d', header="Tc={:.{}f} Eigs={}\n".format(Tc, precision, *eigs))

