import numpy as np
import matplotlib.pyplot as plt
import copy
from scipy import optimize
from scipy import integrate

import sys
sys.path.insert(1, '/home/g/Research/GapFunction/gap/')


def test_temperature_triplet():
    eig = 0.0000000000000001
    Tc_list, Tc2_list = list(), list()
    eigs = np.arange(0,1000,1)
    for i in range(1,1001):
        def func(Tc):
            P = gapFunctions.getWorkMatrix(angular_points,radial_points,Tc,'triplet',0.0)
            P2 = gapFunctions.getWorkMatrix(angular_points,radial_points,Tc,'singlet',0.0)
            w, v = np.linalg.eig(P)
            highEig = gapFunctions.getHighEig(w)
            eig = i / 100 * highEig
            return abs(eig - 1)
        Tc = optimize.least_squares(func, 50, bounds = (0.01, 1000), verbose=0).x[0] 
        #Tc2 = gapFunctions.analyticalApproximate(1/eig,1)
        Tc_list.append(Tc)
        #Tc2_list.append(Tc2)
    plt.plot(Tc_list)
    plt.plot(Tc2_list)
    plt.show()
# Check Function
# Creates eig vs temperature plot and compares to analytical function
def test_temperature(lowEig):
    eigenvalues = [0.01*i*lowEig for i in range(2,1001)]
    invEigenvalues = [1/(0.01*i*lowEig) for i in range(2,1001)]
    numericalValues, analyticalApproximations = list(), list()
    for eig in eigenvalues:
        eig /= lowEig
        p = eig * P
        #gapFunctions.printMatrix(p, 'p')
        temperature,f,eigs,eigenvectors,lowIndex,highIndex, norm = solve(p,N)
        numericalValues.append(temperature)
        Eig = eigs[highIndex]
        analyticalApproximations.append(gapFunctions.analyticalApproximate(1/Eig,N))
        #print(numericalValues[-1][0], eig)
    numericalCol = [x for x in numericalValues]
    analyticalCol = [x[0] for x in analyticalApproximations]
    plt.title('Temperature vs Eigenvalue')
    plt.xlabel('Eigenvalue')
    plt.ylabel('Temperature')
    plt.plot(eigenvalues, numericalCol)
    plt.plot(eigenvalues, analyticalCol)
    plt.show()
    plt.title('Integral vs inverted eigenvalue')
    plt.xlabel('1 / Eigenvalue (Theoretical integral value)')
    plt.ylabel('Integral value')
#    plt.plot(eigenvalues, numericalCol)
    plt.plot(invEigenvalues, numericalCol)
#    plt.plot(eigenvalues, analyticalCol)
    plt.plot(invEigenvalues, analyticalCol)
    plt.show()

# Returns limit where the eigenvalues become so small the code cannot accurately diagonalize the matrix
def test_diagonalization_limit():
    w, v = np.linalg.eig(P)
    eig1 = gapFunctions.pickEigenvalue(w).real
    p = P/10
    w, v = np.linalg.eig(p)
    eig2 = gapFunctions.pickEigenvalue(w).real

    ratio1, ratio2 = 1, 1
    i = 0

    while round(ratio1,10) == round(ratio2,10):
        p /= 10
        w, v = np.linalg.eig(p)
        eig3 = gapFunctions.pickEigenvalue(w).real

        ratio1 = eig2 / eig1
        ratio1 = round(ratio1, 1)
        ratio2 = eig3 / eig2
        ratio2 = round(ratio2, 1)

        eig1 = eig2
        eig2 = eig3

        i += 1
    i -= 1
    return 10**(-i)

# Tests lowest eigenvalue as we expand matrix size
def test_matrix_shape():
    s = list()
    b = list()
    n = list()
    # i is the value n in the matrix shape n*m x n*m
    for i in range(2,83):
        #print(i)
        P = gapFunctions.getWorkMatrix(i, 1, T0, 'default', 0, [0,0,0])
        temperature,f,eigenvalues,eigenvectors,lowIndex, highIndex, norm = solve(P,N)
        #gapFunctions.printMatrix(eigenvectors, 'eigenvectors')
        #for g in eigenvectors.transpose():
        #    plotGapFunction(g)
        lowVal = eigenvalues[lowIndex]*i**2
        highVal = eigenvalues[highIndex]
        s.append(lowVal)
        b.append(highVal)
        n.append(norm)
        #print('Eigenvalue:',val)
        #gapFunctions.printMatrix(P, 'P')
    plt.title("Low Eigenvalue vs shape")
    plt.ylabel("eig")
    plt.plot(range(2,83), s)
    plt.show()
    plt.title("High Eigenvalue vs shape")
    plt.ylabel("eig")
    plt.plot(range(2,83), b)
    plt.show()
    plt.title("Norm vs shape")
    plt.plot(range(2,83), n)
    plt.show()

def test_matrix_shape_2():
    results = list()
    for i in range(2,10):
        V = gapFunctions.getWorkMatrix(i,1,3,'const', 0, [0,0,0])
        gapFunctions.printMatrix(V, 'V')
        eigs, vecs = np.linalg.eig(V)
        eigs, vecs, low, high = gapFunctions.pickEigenvalue(eigs, vecs)
        results.append(eigs[low])
    plt.figure()
    plt.plot(range(2,10), results)
    plt.show()

def test_Q_solve(Tc, size):
    integral = gapFunctions.solveF(Tc)
    row1 = np.array([0.] * size)
    matrix1 = np.ones((size,size))
    for i in range(size):
        E_k = i / size
        if E_k == 0:
            val = 1 / (4*Tc)
        else:
            val = np.tanh(E_k / ((2*Tc))) / (2*E_k)
        val *= 2 / size
        row1[i] = val
        temp_row1 = np.array([row1[i]]*size)
        matrix1[i] = np.matrix.copy(temp_row1)
    gapFunctions.printMatrix(matrix1.transpose(), 'M1')
    print('Integral: ', integral)
    eigs1, vecs1 = np.linalg.eig(matrix1)
    print('List of Eigs: ', eigs1)
    eigs1, vecs1 = gapFunctions.alterVecs(eigs1, vecs1)
    print('High Eig: ', eigs1[0])

def test_g_reduce():
    g = np.array([1,0,0])
    Tc_singlet, w, v, d = solve_complex('singlet',0.0,np.array([0,0,0]))
    temperatures = list()
    print('-=====================================================')
    alpha = 0.0
    Tc, w, v, d = solve_complex('triplet',alpha,g)
    print('Tc_triplet: ', Tc)
    temperatures.append(Tc)
    print("Tc_singlet: ", Tc_singlet)

def temp_alpha():
    t1, t2 = list(), list()
    Tc, eigs, vecs = solve_complex('singlet', 0)
    alphas = list()
    for i in range(100):
        alpha = i / 100
        Tc2, eigs2, vecs2 = solve_complex('triplet', alpha)

        if alpha == 0:
            s = lambda Tc2: 1/(4*Tc2)
            a = 1/(4*Tc2)
        else:
            s = lambda Tc2: 2 / ( 2*Tc2*(np.cosh(alpha / Tc2) + 1) )
            a = 1 * ( 1/(2 * alpha) * np.tanh( alpha / (2*Tc2) ) +  1/(2 * alpha) * np.tanh( alpha / (2*Tc2) ) )

        #Tc = alpha / (2 * np.arctanh(alpha))
        s_list = np.empty([1000])
        for i in range(1,1001):
            s_list[i-1] = s(i / 100)
        if np.amax(s_list) >= 1:
            print('hi')

        if not gapFunctions.check_temp(n,m,'triplet',alpha,Tc2):
            t1.append(Tc)
            print('yo')
            continue
        alphas.append(alpha)
        #print('a: ', a)
        #print('s: ', s)
        #print(eigs2)
        #P = gapFunctions.getWorkMatrix(n,m,Tc,'triplet',alpha)
        ##gapFunctions.printMatrix(P, "P_new")
        #print(Tc,alpha)
        t1.append(Tc)
        t2.append(Tc2)
    plt.figure()
    plt.plot(range(100),t1)
    plt.plot(alphas,t2)
    plt.title('Tc vs alpha')
    plt.xlabel('alpha*100')
    plt.ylabel('Tc')
    plt.show()

def temp_alpha2():
    T_triplet, T_singlet = list(), list()
    alphaMax = 51
    for i in range(alphaMax):
        alpha = i/100
        Tc_triplet, eigs_triplet, vecs_triplet = solve_complex('triplet', alpha)
        Tc_singlet, eigs_singlet, vecs_singlet = solve_complex('singlet', 0)
        a_num, b_num = gapFunctions.Qtransform.a__b(0,Tc_triplet,alpha,np.array([0,1,0]))
        s_num = (b_num-a_num)
        highEig = a_num + (b_num)*(alpha**2-1)
        print(-2*highEig)
        T_triplet.append(Tc_triplet)
        T_singlet.append(Tc_singlet)

    plt.figure(0)
    plt.plot(np.arange(alphaMax)/100,T_triplet)
    plt.plot(np.arange(alphaMax)/100,T_singlet)
    plt.xlabel('alpha')
    plt.ylabel('Superconducting Temperature')
    plt.title('Testing alpha limit')
    plt.show()

def eigPeak(T,alpha):
    if alpha == 0:
        a = 1/(4*T)
        s = 1/(4*T)
    else:
        s = 1 / ( 2*T*(np.cosh(alpha / T) + 1) )
        a = 1/2 * ( 1/(2 * alpha) * np.tanh( alpha / (2*T) ) +  1/(2 * -alpha) * np.tanh( -alpha / (2*T) ) )
    highEig = -a + (-a + s)*(alpha**2 - 1)
    return highEig

def plotEigPeak():
    for j in range(10):
        alpha = j/10
        eigs = list()
        for i in range(1,1001):
            T = i/1000
            a_num, b_num = gapFunctions.Qtransform.a__b(0,T,alpha,np.array([0,1,0]))
            s_num = (b_num-a_num)
            highEig = a_num + (b_num)*(alpha**2-1)
            eig = eigPeak(T,alpha)
            eigs.append(-2*highEig)
        plt.figure()
        plt.plot(np.arange(1000)/1000,eigs)
        plt.xlabel('T')
        plt.ylabel('Eig')
        plt.title('Max Eigenvalue vs temperature for alpha='+str(alpha))
        plt.show()

def solveTc_square_integral(invEig, n):
    from cfg import wc
    def func(Tc):
        def f(x):
            print('x', x)
            E_k = x**2
            if E_k == 0:
                return 1/(4*Tc)
            return np.tanh(E_k/(2*Tc))/(2*E_k)
        sum = 0
        for i in range(n):
            sum += f(i/n*wc) / n * wc
        return abs(sum - invEig)
    return optimize.least_squares(func, 1.0, bounds = (0.01, 1000), xtol = 10**(-14), verbose=0).x[0] 

# Performs integral for value of f(Tc)
# Returns 2 * integral from 0 to wc, works bc the function is symmetric
def solveF(Tc):
    from cfg import wc
    def f(x):
        return np.tanh(x/(2*Tc))/(2*x)
    return integrate.quadrature(f,0,wc)[0] 

def solveTc2(n,m):
    def func(Tc):
        P = getWorkMatrix(n,m,Tc,'singlet',0)
        w, v = np.linalg.eig(P)
        highEig = getHighEig(w)
        print('Singlet Eig: ', highEig)
        #printMatrix(P, 'P')
        return abs(highEig - 1)
    return optimize.least_squares(func, 1.0, bounds = (0.01, 1000), xtol = 10**(-14), verbose=0).x[0] 

def check_temp(n,m,f_type,alpha,Tc):
    P = getWorkMatrix(n,m,Tc,f_type,alpha)
    w, v = np.linalg.eig(P)
    highEig = getHighEig(w)
    print(highEig)
    if round(highEig,8) == 1:
        return True
    return False

# Analytical approximate function given in paper
def analyticalApproximate(invEig,N):
    from cfg import wc
    Tc = 1.134*wc*np.exp(invEig/(-N))
    integral = solveF(Tc)
    #print('Approximate Integral=',solveF(Tc))
    return Tc
    #print('Approximate Critical Temperature=', Tc)
