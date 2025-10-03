# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build and Development Commands

### Building the Project
```bash
# From anywhere if fly-build.sh is added to PATH
./scripts/fly-build.sh        # Standard build (minimal output)
./scripts/fly-build.sh -v     # Verbose build
./scripts/fly-build.sh -vv    # Very verbose build (recommended for first build)
```

The build system uses CMake with Ninja generator and produces:
- **Executable**: `build/bin/fly.x` - Standalone command-line tool
- **Shared Library**: `build/lib/libfly.so` - For Python/Julia bindings

### Running Tests
```bash
fly.x  # Runs all test suites with no arguments
```

Tests must pass before committing. Each category should have a `tests/` folder with `all.cpp` and `all.hpp` following the pattern in `src/objects/tests/all.cpp`.

### Julia Setup
```bash
julia --project=~/Research/FFirefly/jlpkg/Firefly -e 'using Pkg; Pkg.Registry.update(); Pkg.instantiate(); Pkg.precompile()'
```

### Running Calculations
```bash
fly.x < input_file.cfg  # Via config file
# Or use Python wrappers in scripts/ for programmatic access
```

## Architecture Overview

### Multi-Language Integration
Ffirefly is a polyglot codebase designed for computational physics:
- **C**: Main entry point (`src/main.c`)
- **C++**: Core algorithms, objects, and category implementations
- **Python**: Configuration, plotting, high-level wrappers
- **Julia**: Performance-critical many-body calculations (eliashberg, self-energy loops)
- **Fortran**: Legacy integrations (libtetrabz)

The config system (`src/config/`) auto-generates language bindings from `input_variables.py`:
- `write_c.py` → C header/source for global variables
- `write_cpp.py` → C++ config namespace
- `write_py.py` → Python config module
- `write_julia.py` → Julia config module
- `write_f90.py` → Fortran90 config module

### Category-Based Modular Structure

Calculations are organized into **categories** in `src/`:
- `superconductor/` - BCS, Eliashberg theory
- `response/` - Susceptibility, response functions
- `hamiltonian/` - Fermi surface, band structure
- `many_body/` - Vertex, self-energy, renormalization
- `objects/` - Core data structures (CMField, CMData, Bands, Vertex)
- `algorithms/` - Numerical methods, integration

Each category follows the **node wrapper pattern**:
1. Category folder contains implementation files
2. `node.cpp` exports a `{category}_wrapper()` function
3. `node.hpp` declares the wrapper
4. `src/main.c` imports and calls the wrapper based on user's `category` input

Example from `src/superconductor/node.cpp`:
```cpp
extern "C" void superconductor_wrapper() {
    if (calculation == "bcs")
        bcs();
    else if (calculation == "eliashberg")
        eliashberg();
    // ...
}
```

### Sequential Calculations
The `category` field supports `+` separators for chained calculations:
```cfg
category = 'response+vertex+superconductor'
method = 'sparse_ir+FLEX+eliashberg'
```
Each category runs in sequence, allowing output from one to feed into the next.

### Core Data Objects

#### CMData
Structured dataset handler with column-based format:
```
x  y  z  w  n  f
```
- Reads/writes files with headers
- `w` and `n` columns specified during init
- Used throughout for passing physics data between components

#### CMField
Interpolation engine built on CMData:
- Constructs regular mesh grid from `x, y, z` columns
- Binary search interpolation over `w` dimension
- **Automatic periodic boundary conditions** - queries outside grid wrap around
- Base for Field_R, Field_C, Field_VR, Field_VC (real/complex scalar/vector fields)

#### Bands
Provides band structure ε(k):
- Reads from `{prefix}_bands.dat` if available, else uses `.cfg` definition
- **Band indexing starts at 1** (not 0)
- Call: `Bands band; band(n, k)` where `n ∈ [1, nbnd]`

#### Vertex
Two-particle interaction Γ(q, ω):
- Reads from `{prefix}_2PI.dat` or falls back to `.cfg` interaction
- Supports optional labels (spin, valley, orbital)
- Call: `Vertex V; V(q, w, label1, label2)`

### File Naming Convention
All input/output follows: `{prefix}_{filetype}.{ext}`

Examples:
- Bands: `sample_bands.dat` or `sample_bands.h5`
- Susceptibility: `sample_chi.dat`
- Gap function: `sample_gap.dat`
- Vertex: `sample_2PI.dat`

The `prefix` is set in the config file `[CONTROL]` section. Default filetype is now `h5` (HDF5).

## Adding New Features

### Adding a New Category
1. Create folder in `src/{category_name}/`
2. Implement your calculation functions
3. Create `src/{category_name}/node.cpp` and `node.hpp`:
   ```cpp
   extern "C" void {category_name}_wrapper() {
       // Dispatch based on config variables
   }
   ```
4. Add include to `src/main.c`:
   ```c
   #include "{category_name}/node.hpp"
   ```
5. Add dispatcher in main loop (around line 170):
   ```c
   else if (!strcmp(category, "{category_name}"))
       {category_name}();
   ```
6. Define wrapper function before main:
   ```c
   void {category_name}() {
       printf("Starting {Category Name} Calculation\n\n");
       {category_name}_wrapper();
   }
   ```
7. Rebuild with `fly-build.sh`

### Adding Config Variables
Edit `src/config/input_variables.py`:
```python
ALL = {
    "YOUR_CATEGORY": {
        "your_variable": default_value,
    },
}
```
Then run the file to regenerate all language bindings:
```bash
python src/config/input_variables.py
```

Variables become accessible in all languages:
- C: `c_your_variable`
- C++: `your_variable` (in global namespace after `load_cpp_config_wrapper()`)
- Python: `config.your_variable`
- Julia: `your_variable` (after loading config)

### Adding Tests
1. Create `tests/` in your category folder
2. Implement `all.cpp` and `all.hpp` (see `src/objects/tests/all.cpp` as template)
3. Update `num_tests` in `src/main.c` test() function
4. Add your test function call to the `all_tests` array
5. Tests must be fast (<1s) and have verifiable sources (analytical limits or reference papers)

## Important Technical Notes

- **Archive folders are excluded**: CMake skips any `/archive/` directories during compilation
- **OpenMP parallelization**: Uses `num_procs - 1` threads by default
- **ccache enabled**: Speeds up recompilation
- **Python embedded**: Python interpreter initialized in main.c via `start_python()` / `end_python()`
- **Julia called from C++**: Many-body loops use Julia for performance via libjulia
- **Verbosity control**: Set `verbosity = 'high'` in config for detailed output

## Calculation Types and Methods

### Superconductor (`category = 'superconductor'`)
- **calculation**: `bcs`, `eliashberg`, `linearized_eliashberg`
- **method**: `power_iteration`, `diagonalization`, `projection`
- **FS_only**: Compute only on Fermi surface (true) or full BZ (false)
- Output: `{prefix}_gap.dat`

### Response (`category = 'response'`)
- **method**: `sparse_ir` (finite T, dense k-grid), `libtetrabz` (T=0, tetrahedron integration)
- **dynamic**: Include ω-dependence (true) or static ω=0 (false)
- **w_pts**: Number of Matsubara frequencies
- Output: `{prefix}_chi.dat`

### Vertex (`category = 'vertex'`)
- **interaction**: `FLEX` (requires `{prefix}_chi.dat`)
- Output: `{prefix}_vertex.dat`

### DOS (`category = 'DOS'`)
- **method**: `libtetrabz`, `surface_sum`
- Requires dense `k_mesh` for accuracy
- Output: `{prefix}_DOS.dat`

### Many-Body (`category = 'many_body'`)
- Bethe-Salpeter iterations
- **self_consistent**: Iterate to convergence (true) vs single shot (false)
- Implemented in Julia (`src/many_body/many_body_loop.jl`)
