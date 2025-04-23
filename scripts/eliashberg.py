import firefly
import numpy as np
from scipy.interpolate import interp1d

config = {
        'CONTROL': {
            'category': 'superconductor',
            'method': 'eliashberg',
            'prefix': 'sample',
            },
        'SYSTEM': {
            'interaction': 'FLEX',
            'dimension': 2,
            'fermi_energy': -1.2,
            'Temperature': 0.0001,
            'onsite_U': 3.0,
            'nbnd': 1,
            'ibrav': 1,
            },
        'MESH': {
            'k_mesh': [650, 650, 650],
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

def get_slope_diff(x_data, y_data):
    f_interp = interp1d(x_data, y_data, kind='cubic', fill_value='extrapolate')
    x3 = x_data[-1]
    x2 = x_data[-2]
    x1 = x_data[-3]
    y3 = f_interp(x3)
    y2 = f_interp(x2)
    y1 = f_interp(x1)

    m23 = (y3 - y2) / (x3 - x2)
    m12 = (y2 - y1) / (x2 - x1)
    return m23 - m12

def update_dT(T_list, phi_list, dT):
    m = abs(get_slope_diff(T_list, phi_list)) + 1e-4
    r = max(5e-2 / m * dT, 1e-4)
    r = min(r, 5e-2)
    return r

def eliashberg(filename):
    Ti = 1.0e-4
    Tf = 1.0
    dT = 1e-3
    flatline = True
    T_list, phi_list = [], []
    print("Temperature, Max phi")

    for i in range(50):
        config['SYSTEM']['Temperature'] = Ti
        output, error = firefly.launch(config)
        #if(len(error) > 0):
        #    print("Error occured")
        #    print(error)
        #    print("Error lines found: ", len(error))
        match = firefly.grep(output, "Max phi:")
        value = firefly.extract_value(match)

        print(f"{Ti:.5f} {value:.5f}")
        #print(dT)

        T_list.append(Ti)
        phi_list.append(value)

        if i > 2:
            dT = update_dT(T_list, phi_list, dT)
        Ti += dT
        if value < 3e-4 or Ti > Tf:
            break

    data = np.array([T_list, phi_list])
    header = "Temperature Max_phi"
    np.savetxt(filename, data, delimiter=" ", header=header, fmt="%d")

print("Fermi_energy: ", config['SYSTEM']['fermi_energy'])
eliashberg("Phi_v_T.dat")

"""

import ffirefly
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
    output = ffirefly.launch(config)
    match = ffirefly.grep(output, "Max phi:")
    value = ffirefly.extract_value(match)
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
