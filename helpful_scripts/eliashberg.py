import fcode
import numpy as np

config = {
        'CONTROL': {
            'category': 'superconductor',
            'method': 'eliashberg',
            'prefix': 'first',
            },
        'SYSTEM': {
            'interaction': 'FLEX',
            'dimension': 2,
            'fermi_energy': -1.0,
            'Temperature': 0.0001,
            'onsite_U': 4.0,
            'nbnd': 1,
            },
        'MESH': {
            'k_mesh': [100, 100, 100],
            'q_mesh': [10, 10, 10],
            'w_pts': 30000
            },
        'BANDS': {
            'band1': 'tight_binding',
            't0': 1.0
            },
        'SUPERCONDUCTOR': {
            'projections': 's'
            },
        }

def eliashberg(filename):
    Ti = 1e-4
    Tf = 1.0
    phi_max_initial = 0
    flatline = True
    T_list, phi_list = [], []
    print("Temperature, Max phi")

    for i in range(50):
        config['SYSTEM']['Temperature'] = Ti
        output = fcode.launch(config)
        print(output)
        match = fcode.grep(output, "Max phi:")
        value = fcode.extract_value(match)

        if i == 0:
            phi_max_initial = value
        print(round(Ti, 6), round(value, 6))

        T_list.append(Ti)
        phi_list.append(value)
        if phi_max_initial * 0.9 < value:
            Ti += 5e-4
        else:
            Ti += 5e-5
            if flatline:
                Ti -= 5e-4
            flatline = False
        if Ti > Tf or value < phi_max_initial * 0.1:
            break
    data = np.array([T_list, phi_list])
    header = "Temperature Max_phi"
    np.savetxt(filename, data, delimiter=" ", header=header, fmt="%d")

eliashberg("Phi_v_T.dat")

"""

import fcode
import numpy as np
from scipy.optimize import brentq

config = {
        'CONTROL': {
            'category': 'superconductor',
            'method': 'eliashberg',
            'prefix': 'first',
            },
        'SYSTEM': {
            'interaction': 'FLEX',
            'dimension': 2,
            'fermi_energy': -1.0,
            'Temperature': 0.01,
            'onsite_U': 4.0,
            },
        'MESH': {
            'k_mesh': [34, 34, 34],
            'q_mesh': [10, 10, 10],
            'w_pts': 100
            },
        'BANDS': {
            'band': 'tight_binding',
            't0': 1.0
            },
        }

def get_max_phi(T):
    config['SYSTEM']['Temperature'] = T
    output = fcode.launch(config)
    match = fcode.grep(output, "Max phi:")
    value = fcode.extract_value(match)
    return value

T_list, phi_list = [], []

def func(T, max_phi):
    global T_list, phi_list
    value = get_max_phi(T)
    T_list.append(T)
    phi_list.append(value)
    print(T, round(value, 6))
    if value > 2 * max_phi:
        return - max_phi
    return value - max_phi / 10

def eliashberg(filename):
    Ti = 1e-4
    Tf = 1.0
    print("Temperature, Max phi")
    value = get_max_phi(Ti)

    Tc_solution = brentq(lambda T: func(T, value), Tf, Ti, xtol=1e-4)
    print("Tc =", Tc_solution)
    global T_list, phi_list
    data = np.array([T_list, phi_list])
    header = "Temperature Max_phi"
    np.savetxt(filename, data, delimiter=" ", header=header, fmt="%d")

eliashberg("Phi_v_T.dat")
"""
