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
To install the base C/C++ packages, run 

```bash
cd scripts
chmod +x install.sh
./install.sh
```
```
```

This will install the packages sufficient for base functionality. The complete list of packages for all solvers are below.
| Python     | Julia             | C++      | Fortran    | C    |
|:----------:|:-----------------:|:--------:|:----------:|:----:|
| numpy      | PyCall            | g++      | gfortran   | gcc  |
| scipy      | CUDA              | Cmake    | libtetrabz |      |
| matplotlib | FFTW              | BLAS     |            |      |
| h5py       | Roots             | openBLAS |            |      |
| sparse_ir  | SparseIR          | LAPACK   |            |      |
| pandas     | MPI               | LAPACKE  |            |      |
| tbmodels   | PencilFFTs        | Ninja    |            |      |
|            | LoopVectorization | OpenMP   |            |      |
|            |                   | ccache   |            |      |
|            |                   | Boost    |            |      |
|            |                   | pybind   |            |      |
|            |                   | hdf5     |            |      |

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

The easiest way to run Ffirefly is with a `.cfg` input file:

```bash
fly.x < input_file.cfg
```

A basic example is provided in

```bash
sample.cfg
```

Use this file as a reference for the expected input format.

### Input File Structure

Each input file should specify a `category`, which tells Ffirefly what type of calculation to run.

For each category, the input file should also include a matching section:

```cfg
category = CATEGORY_NAME

[CATEGORY_NAME]
...
```

For example, if the category is `superconductor`, then the input file should include a section called

```cfg
[superconductor]
```

This section contains the input variables for that calculation.

### File Prefixes

The `prefix` variable controls the names of files that Ffirefly reads and writes.

Input datasets should follow the format

```bash
prefix_filetype.dat
```

For example, if

```cfg
prefix = hg1201
```

then a density of states file should be named

```bash
hg1201_dos.dat
```

Output files follow the same naming convention.

### Sequential Calculations

Multiple calculations can be run in sequence by joining categories with `+`.

For example:

```cfg
category = bands+dos+superconductor
```

This tells Ffirefly to run each calculation in order.

### Python Wrapper

Ffirefly also provides a Python wrapper for running sequential calculations.

The wrapper includes launchers and data extractors for each calculation type, making it easier to automate workflows and collect results at each step.

This is useful for parameter sweeps, such as calculating a phase diagram over many temperatures and chemical potentials.

An example can be found in

```bash
scripts/eliashberg.py
```

---

## **📚 Developer Guide**

The source code is organized by calculation category.

Each main category has its own folder inside

```bash
src/
```

For example:

```bash
src/superconductor/
src/many_body/
src/hamiltonian/
```

Additional subfolders can be added inside each category as needed.

### Adding a New Category

To add a new calculation category:

1. Create a new folder inside `src/`.

   For example:

   ```bash
   src/new_category/
   ```

2. Add the source files for the new calculation.

3. Add a `node.cpp` file that connects the new category to the main Ffirefly executable.

4. Make sure the new category is called from `main.c`.

5. Recompile the project with

   ```bash
   fly-build.sh
   ```

The `main.c` file reads the input file, determines which category was requested, and calls the appropriate wrapper function.

### Category Nodes

Each category should have a node file that decides which calculation inside that category should be run.

For example, the `superconductor/` folder contains

```bash
src/superconductor/node.cpp
```

This file defines a function called

```cpp
superconductor_wrapper()
```

The `main.c` file calls `superconductor_wrapper()` whenever the input file specifies

```cfg
category = superconductor
```

New categories should follow this same structure.

### Config Variables

New code should interact with the Ffirefly config system.

If your calculation needs new input variables, add them in

```bash
src/config/input_variables.py
```

Your code should read the relevant config variables and modify its behavior accordingly.

Do not hard-code values that should be controlled by the input file.

---

## Testing

Ffirefly is designed to make testing easy across all calculation categories.

Before running large calculations, you should add small tests that check whether the code is working correctly. Good tests can save a lot of time by catching mistakes early.

### Running Tests

To run all existing tests, simply run

```bash
fly.x
```

When `fly.x` is run with no input file, the default test suite is executed.

If any tests fail, Ffirefly will print which tests failed.

### Adding Tests

To add tests for a new category, create a `tests/` folder inside the category folder.

For example:

```bash
src/new_category/tests/
```

Inside this folder, create

```bash
all.cpp
all.hpp
```

Use the format in

```bash
src/objects/tests/all.cpp
```

as a template.

The test file should store the results of each test in a boolean array and call

```cpp
print_test_results
```

to display the results.

Make sure to update

```cpp
num_tests
```

so that it matches the number of tests being run.

Finally, link the test function in `main.c` so that it runs when

```bash
fly.x
```

is executed with no input file.

### Test Requirements

Tests should be quick to run.

A good test should check one or two representative points, not an entire dense mesh.

Each test should also have a clear reference. This can be:

- an analytical limit,
- a known exact result,
- a comparison to a reference paper,
- or a previously validated benchmark.

The goal is to confirm that the calculation is working without making the test suite slow.

---

## Documentation

When adding a new category, update

```bash
User.md
```

The documentation should explain:

- what the category does,
- how the calculation works at a basic level,
- which config variables are required,
- which input files are read,
- which output files are written,
- and how to run a simple example.

New features should not be considered complete until they are documented and tested.
