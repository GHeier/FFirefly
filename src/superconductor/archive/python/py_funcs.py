import numpy as np
from scipy import integrate

t = 1

def epsilon(kx, ky, kz, mu):
    return energy(kx, ky, kz) - mu

def energy(kx, ky=0, kz=0):
    return -2*t*(np.cos(kx) + np.cos(ky) + np.cos(kz))

# e is epsilon, not energy
def f(e, T):
    return 1 / (1 + np.exp(e/T))

def ratio(qx,qy,qz,kx,ky,kz,T,mu):
    e_qk, e_k = epsilon(qx+kx, qy+ky, qz+kz, mu), epsilon(kx,ky,kz, mu)
    if abs(e_qk - e_k) < 0.0001:
        e_qk += 0.001
    return ( f(e_qk, T) - f(e_k, T) ) / ( e_k - e_qk )

def chi(qx, qy, qz, T, mu):
    f = lambda kx, ky, kz: ratio(qx, qy, qz, kx, ky, kz, T, mu) / (2*np.pi)**3
    integral = integrate.tplquad(f, -np.pi, np.pi, -np.pi, np.pi, -np.pi, np.pi, epsabs=0.0005)
    return integral[0]

def test_chi(qx, qy, qz, T, mu):
    sumVal = 0
    N = 55
    k_points = np.linspace(-np.pi, np.pi, N)
    for kx in k_points:
        for ky in k_points:
            for kz in k_points:
                r = ratio(qx, qy, qz, kx, ky, kz, T, mu)
                sumVal += r
    return sumVal / N**3

print(chi(3.14, 3.14, 3.14, 0.25, 0))
print(ratio(3.14, 3.14, 3.14, 1, 0.7, 0.4, 0.25, 0))
print(epsilon(3.14, 3.14, 3.14, 0))
