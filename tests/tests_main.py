import numpy
import matplotlib.pyplot as plt
from alive_progress import alive_bar

import solve_tests

import sys
sys.path.insert(1, '/home/g/Research/GapFunction/gap')
import create_eq
import solve
import cfg


"""
=========================================== Test #1 =================================================

Test analytical approximate to the integration method
"""
def test1():
    cfg.V_name = 'const'
    cfg.wc = 5
    cfg.dim = 2
    cfg.numPts = 4
    cfg.f_type = 'singlet'
    cfg.integration = True
    cfg.coordinates = 'spherical'
    k = create_eq.get_k()
    V_init = create_eq.get_potential(k)
    with alive_bar(983, title='Changing eigenvalue') as bar:
        #x = [i/100 for i in range(18,1001)]
        #exact, approx = list(), list()
        for i in range(18,1001):
            V = (i/100) * V_init
            Tc, eigs, gaps = solve.solve_simple(V)
            Tc2 = solve_tests.analyticalApproximate(1/eigs[0].real,1)
            #exact.append(Tc)
            #approx.append(Tc2)
            if Tc < 0.3*cfg.wc:
                if round(Tc) != round(Tc2):
                    exit('Temperatures do not align with approximation')
            bar()
    #plt.plot(x, exact, label='exact')
    #plt.plot(x, approx, label='approx')
    #plt.legend()
    #plt.show()
    print('Test #1 Passed!')

"""
=========================================== Test #2 =================================================

Test integration method to pure matrix method
"""
def test2():
    cfg.V_name = 'const'
    cfg.dim = 2
    cfg.numPts = 20
    cfg.f_type = 'singlet'
    cfg.coordinates = 'spherical'
    # With integration
    cfg.integration = True
    k = create_eq.get_k()
    V = create_eq.get_potential(k)
    Tc, eigs, gaps = solve.solve_simple(V)

    # Without integration
    cfg.integration = False
    k = create_eq.get_k()
    V = create_eq.get_potential(k)
    Tc2, eigs2, gaps2 = solve.solve(V,k)
    print(Tc,Tc2)


test2()
