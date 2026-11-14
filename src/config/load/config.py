import ast
### Variables ###

#[CONTROL]
category = 'test'
calculation = 'test'
method = 'none'
outdir = './'
debug = False
prefix = 'sample'
verbosity = 'low'
automatic_file_read = True
write_result = True
filetype = 'h5'
mpi_procs = 0

#[SYSTEM]
interaction = 'none'
dimension = 3
celltype = ''
nbnd = 0
fermi_energy = 0.0
num_electrons = 0.0
mu_from_n = False
Temperature = 0.0
cutoff_energy = 0.05
smearing = 0.02
mixing = 0.02
max_iters = 100
qp_weight = 1.0
recurse_level = 0
tol = 0.0001

#[HAMILTONIAN]
hamiltonian = 'tight_binding'
eps_dx2y2 = 0.0
eps_dz2 = 0.0
eps_px = 0.0
eps_py = 0.0
eps_pz = 0.0
delta_dp = 0.0

#[HUBBARD]
U0 = 0.0
U1 = 0.0
J0 = 0.0
J1 = 0.0

#[MESH]
k_mesh = [10, 10, 10]
q_mesh = [10, 10, 10]
w_pts = 100

#[CELL]
cell = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]

#[BRILLOUIN_ZONE]
brillouin_zone = [[6.283185307179586, 0.0, 0.0], [0.0, 6.283185307179586, 0.0], [0.0, 0.0, 6.283185307179586]]

#[BASIS]
states = ['H']
positions = [[0.0, 0.0, 0.0]]

#[BANDS]
band = 'fermi_gas'
eff_mass = 1.0
t0 = 1.0
t1 = 0.0
t2 = 0.0
t3 = 0.0
t4 = 0.0
t5 = 0.0
t6 = 0.0
t7 = 0.0
t8 = 0.0
t9 = 0.0
t10 = 0.0
tz0 = 0.0
tz1 = 0.0
tz2 = 0.0
tz3 = 0.0
tz4 = 0.0

#[SUPERCONDUCTOR]
FS_only = True
num_eigenvalues_to_save = 5
frequency_pts = 0
projections = ''

#[RESPONSE]
dynamic = False

#[MANY_BODY]
self_consistent = False
impurity_solver = 'IPT'
n_cycles = 20000
len_cycles = 200
warmup_cycles = 8000
### End Variables ###
nbnd += 1
### Functions ###

def BZ_from_cell(cell):
    import numpy as np
    a1 = np.array(cell[0])
    a2 = np.array(cell[1])
    a3 = np.array(cell[2])
    b1 = 2*np.pi*np.cross(a2, a3) / np.dot(a1, np.cross(a2, a3))
    b2 = 2*np.pi*np.cross(a3, a1) / np.dot(a2, np.cross(a3, a1))
    b3 = 2*np.pi*np.cross(a1, a2) / np.dot(a3, np.cross(a1, a2))
    return [b1, b2, b3]

def load_config():
    import sys
    import os

    current_file = os.path.abspath(__file__)
    input_file = current_file[:-39] + "build/bin/input.cfg"
    key, value, section = "", "", ""
    got_dimension = False
    index = 0
    with open(input_file, 'r') as f:
        for line in f:
            if "=" not in line and "[" in line and "]" in line:
                section = line.strip()[1:-1]
                index = 0
                continue
            if "=" in line:
                key, value = line.split('=', 1)
            key = key.strip()
            value = value.split('#')[0]  # Remove comments
            value = value.strip()
            value = value.replace("'", "")
            value = value.replace('"', '')
            # Set the variable

#[CONTROL]
            if key == "category":
                global category
                category = value
            if key == "calculation":
                global calculation
                calculation = value
            if key == "method":
                global method
                method = value
            if key == "outdir":
                global outdir
                outdir = value
            if key == "debug":
                global debug
                debug = value == 'true'
            if key == "prefix":
                global prefix
                prefix = value
            if key == "verbosity":
                global verbosity
                verbosity = value
            if key == "automatic_file_read":
                global automatic_file_read
                automatic_file_read = value == 'true'
            if key == "write_result":
                global write_result
                write_result = value == 'true'
            if key == "filetype":
                global filetype
                filetype = value
            if key == "mpi_procs":
                global mpi_procs
                mpi_procs = int(value)

#[SYSTEM]
            if key == "interaction":
                global interaction
                interaction = value
            if key == "dimension":
                global dimension
                dimension = int(value)
                got_dimension = True
            if key == "celltype":
                global celltype
                celltype = value
            if key == "nbnd":
                global nbnd
                nbnd = int(value)
            if key == "fermi_energy":
                global fermi_energy
                fermi_energy = float(value)
            if key == "num_electrons":
                global num_electrons
                num_electrons = float(value)
            if key == "mu_from_n":
                global mu_from_n
                mu_from_n = value == 'true'
            if key == "Temperature":
                global Temperature
                Temperature = float(value)
            if key == "cutoff_energy":
                global cutoff_energy
                cutoff_energy = float(value)
            if key == "smearing":
                global smearing
                smearing = float(value)
            if key == "mixing":
                global mixing
                mixing = float(value)
            if key == "max_iters":
                global max_iters
                max_iters = int(value)
            if key == "qp_weight":
                global qp_weight
                qp_weight = float(value)
            if key == "recurse_level":
                global recurse_level
                recurse_level = int(value)
            if key == "tol":
                global tol
                tol = float(value)

#[HAMILTONIAN]
            if key == "hamiltonian":
                global hamiltonian
                hamiltonian = value
            if key == "eps_dx2y2":
                global eps_dx2y2
                eps_dx2y2 = float(value)
            if key == "eps_dz2":
                global eps_dz2
                eps_dz2 = float(value)
            if key == "eps_px":
                global eps_px
                eps_px = float(value)
            if key == "eps_py":
                global eps_py
                eps_py = float(value)
            if key == "eps_pz":
                global eps_pz
                eps_pz = float(value)
            if key == "delta_dp":
                global delta_dp
                delta_dp = float(value)

#[HUBBARD]
            if key == "U0":
                global U0
                U0 = float(value)
            if key == "U1":
                global U1
                U1 = float(value)
            if key == "J0":
                global J0
                J0 = float(value)
            if key == "J1":
                global J1
                J1 = float(value)

#[MESH]
            if key == "k_mesh":
                global k_mesh
                k_mesh = [int(value.split()[i]) for i in range(3)]
            if key == "q_mesh":
                global q_mesh
                q_mesh = [int(value.split()[i]) for i in range(3)]
            if key == "w_pts":
                global w_pts
                w_pts = int(value)

#[CELL]
            if section == "CELL" and index < 3:
                global cell
                cell.append([float(line.split()[i]) for i in range(3)])
                index += 1

#[BRILLOUIN_ZONE]
            if section == "BRILLOUIN_ZONE" and index < 3:
                global brillouin_zone
                brillouin_zone.append([float(line.split()[i]) for i in range(3)])
                index += 1

#[BASIS]
            if key == "states":
                global states
                states = [value.split()[i] for i in range(3)]
            if section == "BASIS" and index < 3:
                global positions
                positions.append([float(line.split()[i]) for i in range(3)])
                index += 1

#[BANDS]
            if key == "band":
                global band
                band = value
            if key == "eff_mass":
                global eff_mass
                eff_mass = float(value)
            if key == "t0":
                global t0
                t0 = float(value)
            if key == "t1":
                global t1
                t1 = float(value)
            if key == "t2":
                global t2
                t2 = float(value)
            if key == "t3":
                global t3
                t3 = float(value)
            if key == "t4":
                global t4
                t4 = float(value)
            if key == "t5":
                global t5
                t5 = float(value)
            if key == "t6":
                global t6
                t6 = float(value)
            if key == "t7":
                global t7
                t7 = float(value)
            if key == "t8":
                global t8
                t8 = float(value)
            if key == "t9":
                global t9
                t9 = float(value)
            if key == "t10":
                global t10
                t10 = float(value)
            if key == "tz0":
                global tz0
                tz0 = float(value)
            if key == "tz1":
                global tz1
                tz1 = float(value)
            if key == "tz2":
                global tz2
                tz2 = float(value)
            if key == "tz3":
                global tz3
                tz3 = float(value)
            if key == "tz4":
                global tz4
                tz4 = float(value)

#[SUPERCONDUCTOR]
            if key == "FS_only":
                global FS_only
                FS_only = value == 'true'
            if key == "num_eigenvalues_to_save":
                global num_eigenvalues_to_save
                num_eigenvalues_to_save = int(value)
            if key == "frequency_pts":
                global frequency_pts
                frequency_pts = int(value)
            if key == "projections":
                global projections
                projections = value

#[RESPONSE]
            if key == "dynamic":
                global dynamic
                dynamic = value == 'true'

#[MANY_BODY]
            if key == "self_consistent":
                global self_consistent
                self_consistent = value == 'true'
            if key == "impurity_solver":
                global impurity_solver
                impurity_solver = value
            if key == "n_cycles":
                global n_cycles
                n_cycles = int(value)
            if key == "len_cycles":
                global len_cycles
                len_cycles = int(value)
            if key == "warmup_cycles":
                global warmup_cycles
                warmup_cycles = int(value)
            # Finished setting variables
        if not brillouin_zone:
            print("Error: Brillouin zone not specified.")
            brillouin_zone = BZ_from_cell(cell)
        #band = band[1:]

        #nbnd = len(band)
        if outdir[-1] != '/':
            outdir += '/'
    return input_file

def printv(format_string: str, *args):
    if verbosity == "high":
        try:
            print(format_string.format(*args))
        except (IndexError, KeyError) as e:
            print(f"Formatting error: {e}")

# Load the config file on module import
load_config()
