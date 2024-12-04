from scipy.optimize import fmin
from scipy import integrate
import numpy as np
import matplotlib.pyplot as plt
#from alive_progress import alive_bar
import multiprocessing as mp
import tqdm
#from p_tqdm import p_map

import cfg

t = 1
b = 400
N = 40
dk = 2*np.pi/N
#dk = 0.1
#N = int(2*np.pi/dk)
q_points = np.arange(dk, np.pi, dk)
k_points = np.linspace(-np.pi, np.pi, N)

def epsilon(k):
    from cfg import mu
    k = np.append(k,[0, 0])
    return energy(k[0], k[1], k[2]) - mu

def energy(kx, ky=0, kz=0):
    return -2*t*(np.cos(kx) + np.cos(ky) + np.cos(kz))
    return (kx**2 + ky**2 + kz**2)**(1/2)

def k_max():
    from cfg import mu, wc
    #k_max = fmin(lambda kx,ky,kz: abs(epsilon(kx,ky,kz)-wc), [0,0,0], disp=0)
    #return np.pi
    return 1
    return round(np.linalg.norm(k_max),5)

def f_singlet(k, T):
    E_k = epsilon(k)
    if E_k == 0:
        return 1/(4*T)
    return np.tanh(E_k/(2*T))/(2*E_k) 

def g_function(k):
    theta = k[1] + math.pi / 2
    gx, gy, gz = np.cos(theta), np.sin(theta), 0
    g = np.array([gx,gy,gz])
    #g = np.array([1,0,0])
    return g

def f_triplet(k, T):
    from cfg import alpha
    E_k = epsilon(k)
    Q = Qtransform.get_Q_matrix(E_k, T, alpha, g)
    return Q

# e is epsilon, not energy
def f(e):
    if (e==0):
        return 1/2
    elif (e > 0): 
        return 0;
    return 1;
    return 1 / (1 + np.exp((b*(e))))

def ratio(qx,qy,qz,kx,ky,kz):
    e_qk, e_k = epsilon([qx+kx, qy+ky, qz+kz]), epsilon([kx,ky,kz])
    return ( f(e_qk) - f(e_k) ) / ( e_k - e_qk )
    return 1

def chi(qx, qy, qz):
    from cfg import dim
    sumVal = 0
    for kx in (k_points):
        if dim == 1:
            sumVal += ratio(qx, qy, qz,kx,0,0) * dk**1 / (2*np.pi)**1
            continue
        for ky in (k_points):
            if dim == 2:
                sumVal += ratio(qx, qy, qz,kx,ky,0) * dk**2 / (2*np.pi)**2
                continue
            for kz in (k_points):
                print(ratio(qx, qy, qz, kx, ky, kz))
                sumVal += ratio(qx, qy, qz,kx,ky,kz) * dk**3 / (2*np.pi)**3
    return sumVal 

def chi2(qx, qy, qz):
    f = lambda kx, ky, kz: ratio(qx, qy, qz, kx, ky, kz) / (2*np.pi)**3
    integral = integrate.tplquad(f, -np.pi, np.pi, -np.pi, np.pi, -np.pi, np.pi, epsabs=0.01)
    return integral[0]

def get_chi(q):
    q = np.append(q,[0,0])
    q_points = np.loadtxt('/home/g/Research/GapFunction/data/chi/chi_0.1')

    x_index = round((q[0] + np.pi)/dk,5)
    y_index = round((q[1] + np.pi)/dk,5)
    z_index = round((q[2] + np.pi)/dk,5)

    i = int(x_index)*l**2 + int(y_index)*l + int(z_index)
    jx = int(x_index+1)*l**2
    jy = int(y_index+1)*l

    mx = q_points[i+jx][1] - q_points[i][1]
    my = q_points[i+jy][1] - q_points[i][1]
    mz = q_points[i+1][1] - q_points[i][1]

    b = q_points[i][1]

    x = x_index - round(x_index) 
    y = y_index - round(y_index) 
    z = z_index - round(z_index) 

    return [mx*x+b, my*y+b, mz*z+b]

def calculate_chi():
    q_points = np.arange(-np.pi, np.pi, dk)
    l = len(q_points)
    input = np.empty((l**3,4))
    with alive_bar(l**3, title=f'Calculating chi values') as bar:
        for i in range(l):
            for j in range(l):
                for k in range(l):
                    q = np.array([q_points[i], q_points[j], q_points[k]])
                    l = l
                    index = i*l**2 + j*l + k
                    input[index] = np.append(q, chi2(q))
                    bar()
    np.savetxt('/home/g/Research/GapFunction/data/chi/chi_{}'.format(dk), input)

def calculate_chi2():
    q_points = np.linspace(-np.pi, 3.14159, N)
    x, y, z = np.meshgrid(q_points, q_points, q_points)
    x, y, z = x.flatten(), y.flatten(), z.flatten()
    q_points = np.array([[y[i], x[i], z[i]] for i in range(len(x))])

    chi_vals = np.array([p_map(chi2, q_points)])

    #pool = mp.Pool(mp.cpu_count())
    #chi_vals = np.array([pool.map(chi2, q_points)])
    table = np.concatenate((q_points, chi_vals.T), axis=1)
    np.savetxt('/home/g/Research/GapFunction/data/chi/chi_{}'.format(dk), table)

def plot():
    plt.figure(0)
    for mu in [0]:#, -1, -2, -3]:
        cfg.mu = mu
        chi_a_vals,chi_b_vals,chi_c_vals = list(), list(), list()
    #    with alive_bar(len(q_points), title=f'mu = {mu}') as bar:
        for qx in q_points:
            qa = qx*np.array([1,1,1])
            qb = qx*np.array([1,1,0])
            qc = qx*np.array([1,0,0])
            chi_a_vals.append(chi2(qx, qx, qx))
    #            chi_b_vals.append(chi(qb,mu)/(len(k_points)))
    #            chi_c_vals.append(chi2(qx,0,0))
    #            bar()


        plt.plot(q_points, chi_a_vals, label='mu='+str(mu))
#        plt.plot(q_points, chi_b_vals, label='mu='+str(mu))
#        plt.plot(q_points, chi_c_vals, label='mu='+str(mu))
    plt.legend()
    plt.show()

plot()
#print(chi(1,1,1))
