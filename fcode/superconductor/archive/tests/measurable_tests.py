import solve_tests
from alive_progress import alive_it, alive_bar
import sys
import numpy as np

sys.path.insert(1,'/home/g/Research/GapFunction/gap/')
import create_eq
import solve
from settings import *

# Test 0: Compare new matrix to old matrix
# Test 1: Compare k integral to E integral
# Test 2: Compare angular solution to approximate equation
# Test 3: Compare spherical solution to angular solution
# Test 4: Compare cartesian solution to spherical solution
# Test 5: Compare triplet alpha=0 to singlet



# ============================================ Test 0 =============================================
def test0():
    V, k = create_eq.get_potential_spherical()
    P1 = create_eq.getWorkMatrix_ang(numPts,numPts,1,'singlet',0.0,0.0)
    P2 = create_eq.get_matrix(V,k,1)
    if np.array_equal(P1,P2):
        print("Matrices don't match")
        return False
    return True

# ============================================ Test 0 =============================================
def test1():
    with alive_bar(3, title='Test 1') as bar:
        V, k = create_eq.get_potential_spherical()
        bar()
        Tc1, eigs1, gaps1 = solve.solve_simple(V)
        bar()
        Tc2, eigs2, gaps2 = solve.solve_simple_energy(V)
        bar()
    if round(Tc1,6) == round(Tc2,6):
        return True
    else:
        print(Tc1, Tc2)
        print("Test 1 Failed: Temperatures don't match")
    return False

# ============================================ Test 2 =============================================
def test2():
    V, k = create_eq.get_potential()
    bar = alive_it(range(8))
    for _ in bar:
        Tc, eigs, gaps = solve.solve_simple(V)
        Tc2 = solve_tests.analyticalApproximate(1/eigs[0],1)
        if Tc < 0.05:
            if round(Tc,3) != round(Tc2,3):
                exit('Temperatures do not match in approximation')
            text = 'Success!'
        else:
            text = 'Failure'
        print('Tc={}, Tc2={}'.format(Tc, Tc2))
        bar.text(text)
        V /= 2

def test3():
    V, k = create_eq.get_potential()
    Tc1, eigs1, gaps1 = solve.solve_simple(V)
    Tc2, eigs2, gaps2 = solve.solve(V)
    if round(Tc1,3) == round(Tc2,3):
        return True
    return False

test0()
test1()
test2()
