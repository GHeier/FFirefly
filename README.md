# Ffirefly
## **🚀 Welcome**  
Welcome to **The Ffirefly Project**, a Condensed Matter Computational Physics Workspace designed to be easy to both use and extend. Users can run calculations via either a simple input file or intuitive python wrappers. Developers are able to cleanly extend the codebase by adding new models and solvers in whichever language they prefer due to Ffirefly's modular structure.


### Why use Ffirefly?
Let's say an individual wants to run some calculations on the ground state of a material, and lacks an algorithm ready-to-run. They could both learn and code the entire thing from scratch, which is cumbersome. They could also search for a package that performs the relevant computation, but it may be difficult to integrate into their code. Using Ffirefly, all of those issues are alleviated at once. Ffirefly's prewritten algorithms are efficient, can be accessed using multiple languages, and output results in a simple and easy to read manner.

### **🔹 Ffirefly's Current Supported Computational Categories**
- **Bands, Density of States, and Fermi Surface Calculations**
- **Response Functions**
- **Vertex Calculations**
- **Self-Energy Calculations**
- **Superconducting Calculations**

### **🔹 Languages Supported**
- **C/C++** 
- **Python** 
- **Julia** 
- **Fortran** 

#### [User Documentation](./docs/User.md)
#### [Developer Documentation](./docs/Developer.md)

---
### **🔹 Installation**  
The Ffirefly project grants access to a wide variety of extremely powerful algorithms with a simple interface, but depends on many packages that are listed below. Make sure you add /usr/local/lib to your $LD_LIBRARY_PATH.

#### **1️⃣  Required Packages**  
| Python     | Julia         | C++      | Fortran    | C    |
|:----------:|:-------------:|:--------:|:----------:|:----:|
| numpy      | PyCall        | g++      | gfortran   | gcc  |
| scipy      | CUDA          | Cmake    | libtetrabz |      |
| matplotlib | FFTW          | BLAS     |            |      |
| h5py       | Roots         | openBLAS |            |      |
| sparse_ir  | SparseIR      | LAPACK   |            |      |
| pandas     | MPI           | LAPACKE  |            |      |
|            |               | Ninja    |            |      |
|            |               | OpenMP   |            |      |
|            |               | ccache   |            |      |
|            |               | Boost    |            |      |
|            |               | pybind   |            |      |

---
#### **2️⃣ Build Instructions**  
 1) Go to the "scripts" folder and run "./fly-build.sh -vv" to build the code. -v indicates a verbose output, -vv indicates a very verbose output, and a -v option exists for regular verbosity. However, for the first time building, using -vv is recommended in the event of an error. If you are a dev, I recommend setting fly-build.sh to a custom terminal command, so recompilation can be done from outside folders. This may be useful while running tests and material calculations.
 2) To ensure Ffirefly has been properly installed, simply run "fly.x". The default tests will run, and if all pass, then you have downloaded the packages correctly. If not, the package that failed will be listed.

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

### **🔹 Plotting**  
The "plot" folder contains a script, fly-plot.py, which acts as a command-line interface for plotting a given function. This feature is still under development and will allow for easy plotting of results from Ffirefly calculations. It is recommended to make this a custom terminal command so that it can be run from anywhere.

---


## **📚 Developer Guide**
In src/ there are folders for each category of calculation, with subfolders as needed. If you are adding a new category, simply create a new folder in src/ and add a new file for the calculation. The main.c file handles the input file and calls the appropriate calculation function. After your folder has been created, add a "node" that connects main.c to your folder. Make sure to call this node from main.c and compile with fly-build.sh after.
   - An example of this in the superconductor/ folder in the `node.cpp` file. This file has a function called superconductor_wrapper(), which determines the type of calculation to be performed. main.c calls this superconductor_wrapper() function if the category type is "superconductor". New categories should follow this format.
   - If you do add new code to the project, the code must interact with the config file by modifying its behavior based on all relevant variables. Adding config variables is done in `src/config/input_variables.py`.
   - It also must pass tests to confirm that it is working correctly, and to confirm that the code still works upon future development. These tests must have a source, whether it be an analytical limit or a reference paper. This must be quick to run, so simply check that 1 or 2 points are correct rather than calculating an entire dense mesh.

### **🔹 Testing**  
 This project encourages good coding practice by giving easy access to all tests across all projects. By setting up simple tests prior to running large calculations, you can save yourself a lot of time and headache. In Ffirefly, it is easy to set up and run new tests.
 - To run all existing tests, simply run "fly.x". The tests will run, and the ones that fail will be printed. 
 - To add tests to the test suite, make a "tests/" folder in your category folder. In "tests/" create "all.cpp" and "all.hpp", copying the format shown in src/objects/tests/all.cpp. This file stores a bool array of all test results, and runs "print_test_results" to display the results. Replace the tests with your own, and link the function in main.c. Don't forget to update num_tests to the correct number of tests you run.
 - These tests will run when `fly.x` is executed with no arguments.

### **🔹 Documentation**
When adding a new category, describe what it does and how it works in the User.md file. Include the relevant input and outputs, both datafiles and config variables. 

# TODOS
 Add testing for 
    - algorithms
        - integration
        - linear_algebra
    - response
        - susceptibility
    - superconductor
        - bcs
        - eliashberg
        - linearized_eliashberg
    
