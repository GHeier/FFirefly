import fcode
import re

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

def eliashberg():
    Ti = 1e-4
    Tf = 1.0
    print("Temperature, Max phi")
    for i in range(50):
        config['SYSTEM']['Temperature'] = Ti
        output = fcode.launch(config)
        match = fcode.grep(output, "Max phi:")
        value = fcode.extract_value(match)
        print(Ti, round(value, 6))
        Ti *= 2
        if Ti > Tf:
            break

eliashberg()
