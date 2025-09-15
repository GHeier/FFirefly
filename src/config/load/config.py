import ast
### Variables ###

#[CONTROL]
category = 'test'
calculation = 'test'
method = 'none'
outdir = './'
indir = './'
prefix = 'sample'
verbosity = 'low'
automatic_file_read = True
write_result = True
filetype = 'h5'

#[SYSTEM]
interaction = 'none'
dimension = 3
celltype = ''
nbnd = 0
natoms = 0
fermi_energy = 0.0
Temperature = 0.0
onsite_U = 0.0
cutoff_energy = 0.05

#[MESH]
k_mesh = [10, 10, 10]
q_mesh = [10, 10, 10]
w_pts = 100

#[CELL]
cell = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]

#[BRILLOUIN_ZONE]
brillouin_zone = [[6.283185307179586, 0.0, 0.0], [0.0, 6.283185307179586, 0.0], [0.0, 0.0, 6.283185307179586]]

#[ATOMS]
atom = 'X'
position = [0.0, 0.0, 0.0]

#[BANDS]
band = []
band.append('fermi_gas')
eff_mass = []
eff_mass.append(1.0)
t0 = []
t0.append(1.0)
t1 = []
t1.append(0.0)
t2 = []
t2.append(0.0)
t3 = []
t3.append(0.0)
t4 = []
t4.append(0.0)
t5 = []
t5.append(0.0)
t6 = []
t6.append(0.0)
t7 = []
t7.append(0.0)
t8 = []
t8.append(0.0)
t9 = []
t9.append(0.0)
t10 = []
t10.append(0.0)

#[SUPERCONDUCTOR]
FS_only = True
num_eigenvalues_to_save = 1
frequency_pts = 5
projections = ''

#[RESPONSE]
dynamic = False

#[MANY_BODY]
self_consistent = False
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
            value = value.strip()
            value = value.replace("'", "")
            value = value.replace('"', '')
            # Set the variable

#[CONTROL]
            if "category" in key:
                global category
                category = value
            if "calculation" in key:
                global calculation
                calculation = value
            if "method" in key:
                global method
                method = value
            if "outdir" in key:
                global outdir
                outdir = value
            if "indir" in key:
                global indir
                indir = value
            if "prefix" in key:
                global prefix
                prefix = value
            if "verbosity" in key:
                global verbosity
                verbosity = value
            if "automatic_file_read" in key:
                global automatic_file_read
                automatic_file_read = value == 'true'
            if "write_result" in key:
                global write_result
                write_result = value == 'true'
            if "filetype" in key:
                global filetype
                filetype = value

#[SYSTEM]
            if "interaction" in key:
                global interaction
                interaction = value
            if "dimension" in key:
                global dimension
                dimension = int(value)
                got_dimension = True
            if "celltype" in key:
                global celltype
                celltype = value
            if "nbnd" in key:
                global nbnd
                nbnd = int(value)
            if "natoms" in key:
                global natoms
                natoms = int(value)
            if "fermi_energy" in key:
                global fermi_energy
                fermi_energy = float(value)
            if "Temperature" in key:
                global Temperature
                Temperature = float(value)
            if "onsite_U" in key:
                global onsite_U
                onsite_U = float(value)
            if "cutoff_energy" in key:
                global cutoff_energy
                cutoff_energy = float(value)

#[MESH]
            if "k_mesh" in key:
                global k_mesh
                k_mesh = [int(value.split()[i]) for i in range(3)]
            if "q_mesh" in key:
                global q_mesh
                q_mesh = [int(value.split()[i]) for i in range(3)]
            if "w_pts" in key:
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

#[ATOMS]
            if "atom" in key:
                global atom
                atom = value
            if "position" in key:
                global position
                position = [float(value.split()[i]) for i in range(3)]

#[BANDS]
            if "band" in key:
                global band
                band.append(value)
            if "eff_mass" in key:
                global eff_mass
                eff_mass.append(float(value))
            if "t0" in key:
                global t0
                t0.append(float(value))
            if "t1" in key:
                global t1
                t1.append(float(value))
            if "t2" in key:
                global t2
                t2.append(float(value))
            if "t3" in key:
                global t3
                t3.append(float(value))
            if "t4" in key:
                global t4
                t4.append(float(value))
            if "t5" in key:
                global t5
                t5.append(float(value))
            if "t6" in key:
                global t6
                t6.append(float(value))
            if "t7" in key:
                global t7
                t7.append(float(value))
            if "t8" in key:
                global t8
                t8.append(float(value))
            if "t9" in key:
                global t9
                t9.append(float(value))
            if "t10" in key:
                global t10
                t10.append(float(value))

#[SUPERCONDUCTOR]
            if "FS_only" in key:
                global FS_only
                FS_only = value == 'true'
            if "num_eigenvalues_to_save" in key:
                global num_eigenvalues_to_save
                num_eigenvalues_to_save = int(value)
            if "frequency_pts" in key:
                global frequency_pts
                frequency_pts = int(value)
            if "projections" in key:
                global projections
                projections = value

#[RESPONSE]
            if "dynamic" in key:
                global dynamic
                dynamic = value == 'true'

#[MANY_BODY]
            if "self_consistent" in key:
                global self_consistent
                self_consistent = value == 'true'
            # Finished setting variables
        if not brillouin_zone:
            print("Error: Brillouin zone not specified.")
            brillouin_zone = BZ_from_cell(cell)
        #band = band[1:]
        if len(band) != nbnd and nbnd != 1:
            print("Error: Number of bands does not match number of bands specified in input.")
            sys.exit(1)

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
