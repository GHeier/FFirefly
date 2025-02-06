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
            'Temperature': 0.01,
            'onsite_U': 4.0,
            },
        'MESH': {
            'k_mesh': [10, 10, 10],
            'q_mesh': [10, 10, 10],
            'w_pts': 100
            },
        'BANDS': {
            'band': 'tight_binding',
            't0': 1.0
            },
        }

def eliashberg(filename):
    Ti = 1e-4
    Tf = 1.0
    phi_max_initial = 0
    T_list, phi_list = [], []
    print("Temperature, Max phi")

    for i in range(50):
        config['SYSTEM']['Temperature'] = Ti
        output = fcode.launch(config)
        match = fcode.grep(output, "Max phi:")
        value = fcode.extract_value(match)

        if i == 0:
            phi_max_initial = value
        print(Ti, round(value, 6))

        T_list.append(Ti)
        phi_list.append(value)
        if phi_max_initial * 0.9 < value:
            Ti += Tf / 100
        else:
            Ti += Ti / 20
        if Ti > Tf or value < phi_max_initial * 0.1:
            break
    data = np.array([T_list, phi_list])
    header = "Temperature Max_phi"
    np.savetxt(filename, data, delimiter=" ", header=header, fmt="%d")

eliashberg("Phi_v_T.dat")
