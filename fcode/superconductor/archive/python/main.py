from cfg import *
import solve
import create_eq
import potential


"""
===================================== Define File Names =============================================
"""
# Define non-default additions to file names
if f_type == 'singlet':
    extra = 'mu={}'.format(mu)
if f_type == 'triplet':
    extra = 'mu={}_alpha={}'.format(mu,alpha)
extra += '_{}'.format(coordinates)
if integration and coordinates == 'spherical':
    extra += '_simple'
#elif coordinates == 'spherical':
#    extra += '_{}radpts'.format(num_radial_points)

pot_name = "../data/potentials/{}_{:d}D_{:d}pts_{}_pot.dat".format(V_name,dim,numPts,extra)
sol_name = "../data/solutions/{}_{:d}D_{:d}pts_{}_sol.dat".format(V_name,dim,numPts,extra)




"""
===================================== Solve+Save Equation Section ===================================

Block #1 finds potential
Block #2 loads potential from file
Block #3 finds & saves solution
"""
k = create_eq.get_k()
V = create_eq.get_potential(k)
#potential.save(pot_name, V, k)

#V, k = potential.load(pot_name)

Tc, eigs, gaps = solve.solve_simple(V)
print(Tc)
#solve.save(sol_name, Tc, eigs, gaps, k)
