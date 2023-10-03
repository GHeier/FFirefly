import cfg
from cfg import *
import potential
import phys_funcs
from alive_progress import alive_it, alive_bar
from scipy import integrate
from scipy import optimize
from scipy.linalg import eig
import numpy as np
from random import randint
#from prettytable import PrettyTable
import Qtransform

# Takes percentage of radius and angle and converts to actual radial/angle value
# Note: This does not distribute evenly across a surface, as the distance between points increase as the radius value increases
def transform_ang(radius_fraction, angle_fraction):
    wc = 1
    return [ radius_fraction*wc, angle_fraction*2*np.pi ]

def get_f(k, T, mu, f_type, alpha, i, j):
    match f_type:
        case 'singlet':
            return phys_funcs.f_singlet(k, T)
        case 'triplet':
            return phys_funcs.f_triplet(k, T)[i][j]
        case default:
            return 1

# Defines matrix size (m*n x m*n) by P_function
# Divides by number of elements in order to have stable eigenvalues as size increases
def getWorkMatrix_ang(angular_points,radial_points,T,f_type,alpha,mu):
    o = 1
    if f_type == 'triplet':
        o = 2
    scaling_factor = 1 / (angular_points * radial_points * o)
    total_points = angular_points * radial_points * o
    P = np.ones(shape=(total_points,total_points))
    for i in range(radial_points):
        for j in range(angular_points):
            k1 = transform_ang(i / radial_points, j / angular_points)
            for x in range(o):
                for k in range(radial_points):
                    for l in range(angular_points):
                        k2 = transform_ang(k/radial_points, l/angular_points)
                        for y in range(o):
                            match f_type:
                                case 'singlet':
                                    f = phys_funcs.f_singlet(k2,T)
                                case 'triplet':
                                    g = phys_funcs.g_function(k2)
                                    f = phys_funcs.f_triplet(k2,T,alpha,g)[x][y]
                                    #print('g',g)
                                    #print(f*(-1))
                                case default:
                                    f = 1
                            # The 2 is there because we go from 0 to wc instead of -wc to wc
                            # This is allowed because both f(e) functions are symmetric around 0
                            P[i*angular_points*o+j*o+x][k*angular_points*o+l*o+y] = -potential.const(k1, k2) * f * 2 * scaling_factor #* k1[0] * k2[0] * (np.cos(k1[1])*np.cos(k2[1]) + np.sin(k1[1])*np.sin(k2[1]))
    return P

# Converts from 1d matrix index to 3d matrix index
# If x[0] is always 0, then it is 2D
def transform(x_total, points, dimension):
    x = 3*[0]
    x[0] = int(x_total/(points)**(dimension-1))%(points)
    x[1] = int(x_total/(points)**(dimension-2))%(points)
    x[2] = x_total%(points)
    if dimension == 1:
        x[1] = x[2]
    return x

def index_to_k(x):
    from cfg import numPts
    k_max = 1 #0.01*2*np.pi
    k1 = k_max*(2*np.array(x)/(numPts-1)-1)
    if dim < 3:
        k1 = k1[:-1]
    return k1

def index_to_k_spherical(x):
    from cfg import dim, numPts, integration
    k_max = 1
    y = np.array(x) / (numPts-1)
    k = [k_max*y[0], 2*np.pi*y[1], np.pi*y[2]]
    if integration:
        k = k[1:] 
    if dim < 3:
        k = k[:-1]
    return np.array(k)

def spherical_to_cartesian(k):
    from cfg import dim, integration
    edim = dim
    r = k[0] 
    phi = np.pi/2
    if integration:
        edim -= 1
        theta = k[0]
        if dim > 2:
            phi = k[-1]
        r = 1
    elif dim > 1:
        theta = k[1]
    else:
        exit('Coordinate system exchange dimension error')
    new_k = dim*[0]
    new_k[0] = r*np.cos(theta)*np.sin(phi)
    new_k[-1] = r*np.cos(phi)
    new_k[1] = r*np.sin(theta)*np.sin(phi)
    return np.array(new_k)

def get_k():
    from cfg import dim, numPts, integration, coordinates
    eff_dim = dim
    if integration and coordinates == 'spherical':
        eff_dim = dim - 1
    tot_pts = numPts**eff_dim

    k_matrix = list()
    for i in range(tot_pts):
        r1 = transform(i, numPts, eff_dim)
        if coordinates == 'spherical':
            k1 = index_to_k_spherical(r1)
        elif coordinates == 'cartesian':
            k1 = index_to_k(r1)
        else: 
            exit("Wrong coordinates there bub")
        k_matrix.append(k1)

    return k_matrix

def get_V(name, k1, k2, k_matrix):
    match name:
        case 'const':
            return potential.const(k1, k2)
        case 'inv_square':
            return potential.inv_square(k1, k2)
        case 'phonon_coulomb':
            return potential.phonon_coulomb(k1, k2)
        case 'scalapino_singlet':
            return potential.scalapino_singlet(k1, k2, k_matrix)
        case default:
            exit('Unknown Potential')

def get_matrix(V, k, T):
    from cfg import dim, numPts, integration, coordinates
    o = 1
    eff_dim = dim
    if integration and coordinates == 'spherical':
        eff_dim -= 1
    if f_type=='triplet':
        o = eff_dim
    P = np.ones(shape=(o*len(k),o*len(k)))
    for i in range(o*len(k)):
        for j in range(o*len(k)):
            x, y = transform(int(i/o), numPts, eff_dim), transform(int(j/o), numPts, eff_dim)
            k1, k2 = k[int(i/o)], k[int(j/o)]
            if coordinates == 'spherical':
                k1, k2 = spherical_to_cartesian(k1), spherical_to_cartesian(k2)
            f = get_f(k2, T, mu, f_type, alpha, i%o, j%o)
            P[i][j] = f * V[int(i/o)][int(j/o)]
    return P 

def get_potential(k_matrix):
    from cfg import dim, numPts, integration, coordinates
    eDim = dim
    if integration and coordinates == 'spherical':
        eDim -= 1
    V = np.ones(shape=(numPts**eDim,numPts**eDim))
    with alive_bar(V.size, title='Calculating Potential') as bar:
        for i in range(len(k_matrix)):
            k1 = k_matrix[i]
            for j in range(len(k_matrix)):
                k2 = k_matrix[j]
                if coordinates == 'spherical':
                    k1, k2 = spherical_to_cartesian(k1), spherical_to_cartesian(k2)
                V[i][j] = get_V(V_name, k1, k2, k_matrix)
                bar()
    return V
