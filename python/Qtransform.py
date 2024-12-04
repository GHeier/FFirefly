import numpy as np

def a__b(E_k, T, alpha, g):
    mag_g = g[0]**2 + g[1]**2 + g[2]**2
    e_plus = E_k + alpha * mag_g
    e_minus = E_k - alpha * mag_g

    messyAVariable = True
    messySVariable = True

    if E_k == 0:
        if alpha*mag_g == 0:
            a = 1/(4*T)
            s = 1/(4*T)
            messyAVariable = False
            messySVariable = False
        else:
            s = 1 / ( 2*T*(np.cosh(alpha * mag_g / T) + 1) )
            a = 1/2 * ( 1/(2 * e_plus) * np.tanh( e_plus / (2*T) ) +  1/(2 * e_minus) * np.tanh( e_minus / (2*T) ) )
            messyAVariable = False
            messySVariable = False
    if e_minus == 0 and E_k != 0:
        a = 1/2 * ( 1/(2 * e_plus) * np.tanh( e_plus / (2*T) ) +  1/(4*T))
        messyAVariable = False
    if messyAVariable:
        a = 1/2 * ( 1/(2 * e_plus) * np.tanh( e_plus / (2*T) ) +  1/(2 * e_minus) * np.tanh( e_minus / (2*T) ) )
    if messySVariable:
        s = 1 / (2 * E_k) * np.sinh(E_k / T) / ( np.cosh(alpha * mag_g / T) + np.cosh(E_k / T) )

    a*=-1
    b = (a + s)
    #print('E_k, alpha', E_k, alpha)
    #print(a,b,s)
    
    return a, b

def Q_matrix_frtho(a,b,g):
    gx, gy, gz = g / (g[0]**2 + g[1]**2 + g[2]**2)
    q11 = a+b*(gx**2-1)
    q12 = b*gx*gy
    q13 = b*gx*gz
    q21 = b*gx*gy
    q22 = a+b*(gy**2-1)
    q23 = b*gy*gz
    q31 = b*gx*gz
    q32 = b*gy*gz
    q33 = a+b*(gz**2-1)

    Q = np.array([ [q11,q12,q13], [q21,q22,q23], [q31,q32,q33] ])
    Q = np.array([ [q11,q12], [q21,q22] ])
    return Q

def shapeQ(Q):
    newQ = np.eye(3)
    for i in range(len(Q)):
        for j in range(len(Q[i])):
            newQ[i][j] = Q[i][j]
    return newQ

def get_Q_matrix(E_k, T, alpha, g):
    a,b = a__b(E_k, T, alpha, g)
    Q = Q_matrix_frtho(a,b,g)
    #print(g)
    #print(Q*-2/5)
    return Q

