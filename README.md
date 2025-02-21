# **QuasiCode**
---

## **🚀 Welcome**  
Welcome to **Quasi**, a computational physics toolkit designed to be easy to both use and extend. Users interact via a simple input file and easy python wrappers. Developers can easily extend the codebase by adding new models and solvers in whatever language they prefer.

### **🔹 Languages Supported**
- **Fortran** 
- **C** 
- **C++** 
- **Python** 
- **Julia** 

### **🔹 Current Calculation Categories**
- Superconducting Gap calculations, both BCS and Eliashberg
-- cateogry name: superconductor
-- methods: bcs, eliashberg
-- calculation: `diagonalization` of exact gap functions or `projection` of gap functions onto a basis set
-- FS_only: `true` or `false` to calculate the gap only on the Fermi surface or the entire Brillouin zone
- Compute **self-energy (Σ)** and **2PI vertices (Γ) (with possible corrections)**  
- Response functions (polarization, etc)
-- category name: response
-- methods: libtetrabz(analytic tetrahedra integration), sparse_ir (sparse matsubaras with convolution theorem)
-- wpts: number of Matsubara frequencies to use. wpts=1 calculates iv=0
- Band structure/density of states/Fermi surface calculations and displays

---

## **📖 User Guide**  
The simplest way to run quasi.x is to create an input file and run it with the command `quasi.x < input_file.cfg`. An example can be seen in the file `sample.cfg`. It is recommended to look at this file to understand the syntax of the input file. It is inspired from the Quantum Espresso and LmtART input file formats.
 - `category` indicates the type of calculation to be performed, with a section [CATEGORY_NAME] in the file for inputs specific to that category. Prefix is the filename prefix for files read and written by this program. Files are saved in the format 'prefix_filetype.dat', so for instance a density of states calculation would save to 'prefix_dos.dat'.


### **🔹 Installation**  
The upside of a package like this is that it grants access to a wide variety of extremely powerful algorithms with a very simple interface. The downside is that because it contains so many packages, it can be a bit of a pain to get everything set up. The packages required are listed below.

#### **1️⃣  Required Packages**  
| Python     | Julia         | C++    | Fortran    | C    |
|:----------:|:-------------:|:------:|:----------:|:----:|
| numpy      | PyCall        | g++    | gfortran   | gcc  |
| scipy      | CUDA          | Ninja  | libtetrabz |      |
| matplotlib | FFTW          | BLAS   |            |      |
| h5py       | Roots         | LAPACK |            |      |
| sparse_ir  | SparseIR      | CMake  |            |      |
| pandas     | LinearAlgebra | OpenMP |            |      |
|            | Printf        | ccache |            |      |


#### **2️⃣ Build Instructions**  
 - Once the above have been downloaded, simply go to the "helpful_scripts" folder and run "./qbuild.sh -vv" to build the code. -v indicates a verbose output, -vv indicates a very verbose output. I recommend setting qbuild.sh to a terminal command, so that you can recompile from outside folders (useful when running tests and material calculations).
 - To run all tests, simply run "quasi.x". The tests will run, and if all pass then you have downloaded the packages correctly. If not, the package that failed will be printed to the screen.

## **📚 Developer Guide**
 - This code is designed to encourage good coding practice but to not impede the developer. In src/ there are folders for each category of calculation, with subfolders as needed. If you are adding a new category, simply create a new folder in src/ and add a new file for the calculation. The main.c file handles the input file and calls the appropriate calculation function. After your folder has been created, add a "node" that connects main.c to your folder. Make sure to call this node from main.c and compile with qbuild.sh after.
 -- An example of this in the superconductor/ folder in the `node.cpp` file. This file has a function called superconductor_wrapper(), which determines the type of calculation to be performed. main.c calls this superconductor_wrapper() function if the category type is "superconductor". New categories should follow this format.

### **🔹 Testing**  
 By far, the most time spent coding is actually spent debugging. This project encourages good coding practice by giving easy access to all tests across all projects. By setting up simple tests prior to running large calculations, you can save yourself a lot of time and headache. It is easy to set up, and a good habit to get into.
 - To run all tests, simply run "quasi.x". The tests will run, and the ones that fail will be printed. 
 - To add tests to the test suite, make a "tests/" folder in your category folder. In "tests/" create "all.cpp" and "all.hpp". A good example is in src/objects/tests/all.cpp, and should be copied to your cateogry folder. This file stores a bool array of all test results, and runs "print_test_results" to display the results. Simply copy this file, replace the tests with your own, and link the function in main.c. Don't forget to update num_tests to the correct number of tests you run.

### **🔹 Documentation**
Documentation is important, but it is also important to not let it get in the way of coding. If adding a new category, describe what it does and how it works in the README.md file. Include the relevant input and outputs, both datafiles and config variables. 

### ** 🔹 Plotting**  
 - Plotting is handled in the "qplot" folder. This folder contains a "base" python script that acts as a command-line interface to the plotting functions. Still under development, this allows for easy plotting of results. Just as qbuild.sh, I recommend making this a command line argument so that it can be run from anywhere.

