# Ffirefly
## **🚀 Welcome**  
Firefly is a computational physics workspace for condensed matter research. It provides ready-to-use algorithms, simple input-driven workflows, and a modular structure that makes it easy for anyone to add new models and solvers, connect them to existing calculations, and build large computational pipelines.

### **🔹 Current Computational Areas**
- **Electronic Structure**: Bands, Density of States, and Fermi Surface Calculations
- **Many Body Physics**: Greens functions, self-energies, vertices and response functions
- **Superconductivity**: BCS and Eliashberg solvers

### **🔹 Development Languages Supported**
- **C/C++** 
- **Python** 
- **Julia** 
- **Fortran** 

#### [User Documentation](./docs/User.md)
#### [Developer Documentation](./docs/Developer.md)

---
### **🔹 Installation**  
The Ffirefly project consists of the base executable and additional methods. The base package uses C/C++, with the python and julia packages used for various methods.

Make sure you add /usr/local/lib to your $LD_LIBRARY_PATH for ease of compilation. 

#### **1️⃣  Required Packages**  
For the base install to work, the below C/C++ packages are needed. Install before following the Build Instructions
| C/C++      |
|:--------:|
| gcc      |
| g++      |
| Cmake    |
| BLAS     |
| openBLAS |
| LAPACK   |
| LAPACKE  |
| Ninja    |
| OpenMP   |
| ccache   |
| Boost    |
| pybind   |
| hdf5     |

The complete list of packages required for the various methods are below. These are not needed for base functionality. Install as needed, after confirming download works.

| Python     | Julia             | Fortran    |
|:----------:|:-----------------:|:----------:|
| numpy      | PyCall            | gfortran   |
| scipy      | CUDA              | libtetrabz |
| matplotlib | FFTW              |            |
| h5py       | Roots             |            |
| sparse_ir  | SparseIR          |            |
| pandas     | MPI               |            |
| tbmodels   | PencilFFTs        |            |
|            | LoopVectorization |            |
|            |                   |            |

---

#### 2. Build Instructions

1. From the `scripts/` directory, build Ffirefly with

   ```bash
   ./fly-build.sh -vv
   ```

2. After building, confirm all tests pass by running

   ```bash
   cd build/bin
   ./fly.x
   ```

3. Some Ffirefly modules use Julia. To precompile the Julia environment, run

   ```bash
   julia --project=/path/to/FFirefly/jlpkg/Firefly -e 'using Pkg; Pkg.Registry.update(); Pkg.instantiate(); Pkg.precompile()'
   ```

   Replace `/path/to/FFirefly` with the path to your local Ffirefly installation.

> **Developer notes:** You should make fly.x custom terminal command. In bash, this can be done with 

   ```bash
   nano ~/.bashrc
   alias fly.x="/pathto/fly.x"
   ```

   Replace /pathto with the path to fly.x in build/bin

> If you rebuild often, it may be useful to add `fly-build.sh` as a custom terminal command as well. This lets you recompile Ffirefly from any directory while testing code or running material calculations.

> For fly-build.sh, the `-v` and `-vv` flags control the verbosity of the build output. For a first build, `-vv` is recommended because it makes compilation errors easier to diagnose.


---

## **📖 User Guide**  
The simplest way to use Ffirefly is to create a .cfg input file and run fly.x with the command `fly.x < input_file.cfg`. An example input can be seen in the file `sample.cfg`, which shows the proper formatting. 

    `category` indicates the type of calculation to be performed, 
    with a separate section called [CATEGORY_NAME] for inputs.
    "prefix" is the filename prefix for files read and written by this program. 
    Any custom datasets to be read by Ffirefly MUST be named in the format of 'prefix_filetype.dat', so a density of states calculation should be called 'prefix_dos.dat' and so on. Output files are saved following the same naming convention.
    Sequential calculations may be run using '+' signs between each specified category.

A python wrapper is also available to run multiple sequential calculations.  
This wrapper comes with a launcher and data extracters for every calculation type, making it easy to create sets of calculations and extract the data at every step. This can be useful, for instance, to calculate the phase across a range of temperatures and chemical potentials. An example of this can be seen in scripts/eliashberg.py.

---

## **📚 Developer Guide**
In src/, there are folders for each category of calculation, with subfolders as needed. If you are adding a new category, simply create a new folder in src/ and add a new file for the calculation. After your folder has been created, add a "node" that connects main.c to your calculation folder and make sure to call it within main.c. Finally, recompile the project with fly-build.sh. The main.c file handles the input file and calls the appropriate calculation function.
An example of this in the superconductor/ folder in the `node.cpp` file. This file has a function called superconductor_wrapper(), which determines the type of calculation to be performed. main.c calls this superconductor_wrapper() function if the category type is "superconductor". New categories should follow this format.
If you do add new code to the project, the code must interact with the config file by modifying its behavior based on all relevant variables. Adding config variables is done in `src/config/input_variables.py`.
   - It also must pass tests to confirm that it is working correctly and that the code still works upon future development. These tests must have a source, whether it be an analytical limit or a reference paper. This must be quick to run, so simply check that 1 or 2 points are correct rather than calculating an entire dense mesh.

### **🔹 Testing**  
 This project encourages good coding practice by giving easy access to all tests across all projects. By setting up simple tests prior to running large calculations, you can save yourself a lot of time and headache. In Ffirefly, it is easy to set up and run new tests.
 - To run all existing tests, simply run "fly.x". The tests will run, and the ones that fail will be printed. 
 - To add tests to the test suite, make a "tests/" folder in your category folder. In this new folder, create "all.cpp" and "all.hpp", copying the format shown in src/objects/tests/all.cpp. This file stores a bool array of all test results, and runs "print_test_results" to display them. Replace the tests with your own, and link the function in main.c. Don't forget to update num_tests to the correct number of tests you run.
 - These tests will run when `fly.x` is executed with no arguments.

### **🔹 Documentation**
When adding a new category, describe what it does and how it works in the User.md file. Include the relevant input and outputs, both datafiles and config variables. 
