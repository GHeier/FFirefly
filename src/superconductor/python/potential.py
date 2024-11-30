from phys_funcs import chi
import numpy as np

def const(k1, k2):
    return 1

def inv_square(k1, k2):
    k = k1 - k2 
    return 1 / (1 + np.linalg.norm(k))

def phonon_coulomb(k1, k2):
    wc = 1
    wd = 1
    wp = 1
    # Define k
    k1x, k1y = k1[0]*np.cos(k1[1]), k1[0]*np.sin(k1[1])
    k2x, k2y = k2[0]*np.cos(k2[1]), k2[0]*np.sin(k2[1])
    k1z, k2z = 0, 0

    # Phonon potential
    diff = ( (k1x-k2x)**2 + (k1y-k2y)**2 )**(1/2)
    Vp = 1/3
    if diff != 0:
        Vp = 1/(1 + 2 * (k1x-k2x)**2 / diff**2 ) 
    if k1[0]**2 > wd or k2[0]**2 > wd:
        Vp = 0

    # Coulomb Potential
    Vc = 1 / (1 + diff**2)
    if k1[0]**2 > wp or k2[0]**2 > wp:
        Vc = 0

    return Vp + Vc

def scalapino_singlet(k1, k2, k_matrix):
    t = 1
    U = 1
    # Define k
    k_add = k1+k2
    k_sub = k1-k2
    Vs = U**2 * chi(k_add, k_matrix) / ( 1 - U*chi(k_add, k_matrix) ) + U**3 * (chi(k_sub, k_matrix))**2 / (1 - U**2 * (chi(k_sub, k_matrix))**2)
    
    return Vs

def f_singlet(k,T):
    E_k = k[0]#**2
    if E_k == 0:
        return -1/(4*T)
    return -np.tanh(E_k/(2*T))/(2*E_k)

def f_triplet(k, T, alpha, g):
    E_k = k[0]**2
    #print(type(E_k), type(T), type(alpha), type(g))
    Q = Qtransform.get_Q_matrix(E_k, T, alpha, g)
#    print('g',g)
#    print('Q',Q)
    return Q

def g_function(k):
    theta = k[1] + math.pi / 2
    gx, gy, gz = np.cos(theta), np.sin(theta), 0
    g = np.array([gx,gy,gz])
    #g = np.array([1,0,0])
    return g

def save(name, V, k):
    new_k = np.empty((len(k)**2,len(k[0])*2))
    for i in range(len(k)):
        for j in range(len(k)):
            new_k[i*len(k)+j] = np.concatenate([k[i],k[j]])
    new_V = V.flatten()
    new_V = new_V[..., None]
    table = np.concatenate([new_k, new_V], axis=1)
    np.savetxt(name, table)

def load(name):
    table = np.loadtxt(name)
    V = table[:,-1]
    k = table[:,:-1]

    L = int(np.sqrt(len(V)))
    V = np.reshape(V, (L,L))

    eff_dim = int(len(k[0])/2)
    k = k[:numPts**eff_dim,eff_dim:]
    return V, k
