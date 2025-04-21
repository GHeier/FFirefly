# Ffirefly
## **🚀 Welcome**  
Welcome to **The Ffirefly Project**, a Condensed Matter Computational Physics Workspace designed to be easy to both use and extend. Users interact via a simple input file and easy python wrappers. Developers can easily extend the codebase by adding new models and solvers in whatever language they prefer.

### Why use Ffirefly?
Scenario: Professor X wants to run some calculations to determine the ground state of a material, and lacks an algorithm ready-to-run. Professor X has 3 options:
1) Learn the algorithm/theory, code it from scratch, and then integrate it into their code (really slow)
2) Find a package that performs this computation, become familiar with it, and then integrate it into their code (slow)
3) Use Ffirefly, with prewritten calculations, and near-automatic code integration (fast)

### **🔹 Ffirefly's Current Supported Computational Categories**
- **Bands, Density of States, and Fermi Surface Calculations**
- **Response Functions**
- **Vertex Calculations**
- **Self-Energy Calculations**
- **Superconducting Calculations**

#### [User Documentation](./docs/User.md)
#### [Developer Documentation](./docs/Developer.md)

### **🔹 Languages Supported**
- **C/C++** 
- **Python** 
- **Julia** 
- **Fortran** 

---

## **📖 User Guide**  
The simplest way to use Ffirefly is to create a .cfg input file and run fly.x with the command `fly.x < input_file.cfg`. An example input can be seen in the file `sample.cfg`. Look at this file to understand the proper syntax. Refer to the [User Documentation](./docs/User.md) for a more in depth explanation.

    `category` indicates the type of calculation to be performed, 
    with a section [CATEGORY_NAME] in the file for inputs specific to that category. 
    Prefix is the filename prefix for files read and written by this program. 
    Files are saved in the format 'prefix_filetype.dat', so for instance a 
    density of states calculation would save to 'prefix_dos.dat'.

A python wrapper is available to run multiple sequential calculations. This may also be accomplished in the main executable by the use of '+' signs in between sequential calculation categories in the `category` input line. 
The python wrapper comes with launchers and data extracters for every calculation type, so it is easy create sets of calculations and extract the data at every step. This can be useful, for instance, to calculate the phase across a range of temperatures and chemical potentials. An example of this can be seen in scripts/eliashberg.py.

---

### **🔹 Plotting**  
Plotting is handled in the "plot" folder. This folder contains a "fly-plot" python script that acts as a command-line interface to the plotting functions. Still under development, this allows for easy plotting of results. Just as fly-build.sh, it is recommended to make this a custom terminal command so that it can be run from anywhere.

---

### **🔹 Installation**  
This package grants access to a wide variety of extremely powerful algorithms with a very simple interface. This project depends on many packages to run properly, which are listed below. Make sure you add /usr/local/lib to your $LD_LIBRARY_PATH.

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





#### **2️⃣ Build Instructions**  
 - Once the above have been downloaded, simply go to the "scripts" folder and run "./fly-build.sh -vv" to build the code. -v indicates a verbose output, -vv indicates a very verbose output. The first time building, using -vv can be helpful in the event of an error. I recommend setting fly-build.sh to a custom terminal command, so recompilation can be done from outside folders, which can be useful when running tests and material calculations.
 - To ensure Ffirefly has been properly installed, simply run "fly.x". The tests will run, and if all pass then you have downloaded the packages correctly. If not, the package that failed will be listed.
 - Adding tests and re-running this is also an excellent way to test your code while creating projects

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

### **🔹 Documentation**
When adding a new category, describe what it does and how it works in the User.md file. Include the relevant input and outputs, both datafiles and config variables. 

