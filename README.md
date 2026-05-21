# Ffirefly

- [🚀 Welcome](#-welcome)
- [📦 Installation](#-installation)
- [📖 User Guide](#-user-guide)
- [📚 Developer Guide](#-developer-guide)
- [🧪 Testing](#-testing)
- [📝 Documentation](#-documentation)


---

## **🚀 Welcome**  
Firefly is a computational physics workspace for condensed matter research. It provides ready-to-use algorithms, simple input-driven workflows, and a modular structure.

This structure makes it easy for anyone to add new models and solvers, connect them to existing calculations, and build large computational pipelines.

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
## **📦 Installation**  
The Ffirefly project consists of the base executable and additional methods. The base package uses C/C++, with the python and julia packages used for various methods.

Make sure you add /usr/local/lib to your $LD_LIBRARY_PATH for ease of compilation. 

#### 1. Required Packages
To install the base C/C++ packages, run 

```bash
cd scripts
chmod +x install.sh
./install.sh
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
| triqs      | LoopVectorization | OpenMP   |            |      |
| triqs_tprf |                   | ccache   |            |      |
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
   cd ../build/bin
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

A basic example is provided in `sample.cfg`. Use this file as a reference for the expected input format.

### Input File Structure

Each input file should specify a `category`, `calculation` type and `method` choice. This tells the executable what project to run.

### File Prefixes

The `outdir` and `prefix` variables control the names of files that Ffirefly reads and writes. Input & output datasets follow the format

```bash
outdir/prefix_filetype.h5
```

For example, if

```cfg
prefix = hg1201, outdir='./data'
```

then the density of states file is named

```bash
data/hg1201_dos.h5
```
This convention must be followed when saving results so that Firefly projects know where to look for input data. Functions for aving and reading data are supplied, with a breakdown in the [Developer Documentation](./docs/Developer.md).


### Sequential Calculations

Multiple calculations can be run in sequence by joining categories with `+`.

For example:

```cfg
category = bands+dos+superconductor
```

This tells Ffirefly to run each calculation in order.

### Python Wrapper

Ffirefly also provides a Python wrapper for running sequential calculations. The wrapper includes launchers and grep functions to read command line output. 

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

### Adding a New Category

To add a new calculation category:

1. Go to `src/config` and edit `categories.py`, which stores all the categories, calculations, and methods.

2. Add the names of your category, calculation, and/or method.

3. Run `python categories.py` to add your project to Firefly. Folders and files will be created for you.

4. Run `fly-build.sh` to recompile base package

5. Navigate to `src/{category}/{calculation}{method}` to see the setup of your project folder.

### Category Nodes

Each project folder has `run`, `README.md`, and `tests/`. The `run` file is the one that is executed when calling `fly.x`, the `README.md` is your documentation, and `tests/` contains the testing file. 

`tests/` is helpful for personal testing of your code, and for others to confirm that your code works properly on their machine.

> Developer note: Run the `test` file in `tests/`, it's really helpful!

### Config Variables

If your calculation needs new input variables, add them in

```bash
src/config/input_variables.py
```

and then run `python input_variables.py`

Your code has access to config variables (from the input file). There are examples of how to access these config variables from your code in every newly created project.

---

## **🧪 Testing**

Ffirefly is designed to make testing easy across all calculation categories. Good tests can save a lot of time by catching mistakes early.

### Running Tests

When `fly.x` is run with no input file, the default test suite is executed. If any tests fail, Ffirefly will print which tests failed.

To test beyond the default suite, navigate to a project's test folder and run the test file. It will output success or failure based on the test conditions set by the developer.

> Developer note: This is very helpful for debugging, being able to have a constant set of tests to check when you change your code.

### Test Requirements

Ideal tests are quick to run. A good test should check one or two representative points, not an entire dense mesh.

Each test should also have a clear reference. This can be:

- an analytical limit,
- a known exact result,
- a comparison to a reference paper,
- or a previously validated benchmark.

The goal is to confirm that the calculation is working without making the test suite slow.

---

## **📝 Documentation**

When adding a new category, update your local `README.md` file. The `Quick Description` will be seen in the `User.md` file guide.

The rest of the `README.md` is seen upon viewing, and should explain:

- what the code does,
- how the calculation works at a basic level,
- which config variables are required,
- which input files are read,
- which output files are written,
- and any dependencies (ie numpy)

If all tests meet the conditions above, and documentation is filled out, your code may be added to the main repository for all to use.
