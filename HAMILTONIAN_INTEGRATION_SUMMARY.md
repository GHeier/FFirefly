# Hamiltonian Object Integration Summary

## Completed Tasks

### 1. Preloaded Hamiltonians Library ✓

Created a comprehensive preloaded Hamiltonians library in `src/hamiltonian/`:

**Files Created:**
- `preloaded_hamiltonians.hpp` - Header with declarations
- `preloaded_hamiltonians.cpp` - Implementation of all models
- `README_PRELOADED.md` - Complete documentation

**Available Models:**
1. **Simple Cubic (SC)** - Single-band tight-binding with up to 3rd nearest-neighbor hopping
2. **Face-Centered Cubic (FCC)** - 12 nearest neighbors at face centers
3. **Body-Centered Cubic (BCC)** - 8 nearest neighbors at corners
4. **Fermi Gas** - Free electron gas with quadratic dispersion
5. **Emery Model** - 3×3 multi-band model for cuprate superconductors

All models return `vector<vector<complex<double>>>` format for consistency.

### 2. Hamiltonian Class Integration ✓

Modified `src/objects/CMField/hamiltonian.cpp` and `hamiltonian.hpp` to:
- Automatically use preloaded models when HDF5 file not found
- Use the `method` config variable to select lattice type (sc/fcc/bcc/fermi_gas/emery)
- Add `use_preloaded` and `file_found` boolean flags
- Implement `evaluate_preloaded()` method for model evaluation

**Integration Details:**
- When `automatic_file_read = true` and file not found, checks `method` config
- If `method` matches a known preloaded model, uses that instead of failing
- Seamlessly converts between complex<double> (preloaded) and complex<float> (Field_CM)
- Parameters extracted from config: t1, t2, t3, fermi_energy

### 3. C++ Tests ✓

Created comprehensive test suite in `src/hamiltonian/tests/`:

**Files Created:**
- `test_preloaded.cpp` - Individual test functions for all models
- Modified `all.cpp` - Integrated 12 new tests into FFirefly test framework

**Tests Implemented:**
- SC model at Gamma, X, and R points
- BCC model at Gamma point
- Fermi gas dispersion
- Emery model Hermiticity and dimensions
- Lattice type string parsing
- Hamiltonian class with different models
- Matrix element access and validation

**Test Results:** All 14 Hamiltonian tests passing ✓

### 4. C++ Exports ✓

Added exports in `src/module/exports/cpp_export.cpp`:

**New Export Functions:**
- `Hamiltonian_operator_export_list()` - Batch evaluation for multiple k-points
- `Hamiltonian_uses_preloaded()` - Check if using preloaded model
- `Hamiltonian_file_found()` - Check if loaded from file

**Note:** `Hamiltonian_export0()` and `Hamiltonian_operator_export0()` already existed.

### 5. Python Bindings ✓

Updated `src/module/imports/cpp_imports.py`:

**Enhancements:**
- Added list k-point evaluation support
- Added `uses_preloaded` property
- Added `file_found` property
- Supports both single k-point and list of k-points
- Returns numpy arrays of complex64

**Python Test Results:**
```
=== Testing SC Preloaded Hamiltonian ===
Uses preloaded: True
File found: False
H(Gamma) = -7.400000+0.000000j
✓ Gamma point energy correct
✓ List evaluation works
✓ SC Hamiltonian test passed!

=== Testing Emery Preloaded Hamiltonian ===
H shape: (3, 3)
✓ Emery model is 3x3 and Hermitian
✓ Emery Hamiltonian test passed!

=== Testing Multiple Models ===
SC    at Gamma: E = -10.099999
BCC   at Gamma: E = -12.500000
FCC   at Gamma: E = -10.500000
✓ All models work
```

### 6. Julia Bindings ✓

Updated `src/module/imports/cpp_imports.jl`:

**Enhancements:**
- Added `uses_preloaded()` function
- Added `file_found()` function
- Added list k-point evaluation method
- Supports both single k-point and vector of k-points
- Returns Vector{Matrix{ComplexF32}}

**Exported Functions:**
- `Hamiltonian` struct
- `uses_preloaded(ham)`
- `file_found(ham)`

## Usage Examples

### C++
```cpp
#include "objects/CMField/hamiltonian.hpp"
#include "config/load/cpp_config.hpp"

// Set config
method = "sc";
t1 = 1.0;
t2 = 0.2;
fermi_energy = 0.0;

// Create Hamiltonian
Hamiltonian ham;

// Check status
if (ham.use_preloaded) {
    cout << "Using preloaded model: " << method << endl;
}

// Evaluate at k-point
Vec k;
k.x = M_PI/2; k.y = M_PI/2; k.z = 0.0;
auto H = ham(k);  // Returns vector<vector<complex<float>>>
```

### Python
```python
from cpp_imports import Hamiltonian, load_config
import numpy as np

# Configure for SC model
load_config("input.cfg")  # Sets method = 'sc', t1, etc.

# Create Hamiltonian
ham = Hamiltonian()

# Check status
print(f"Uses preloaded: {ham.uses_preloaded}")

# Single k-point
k = [np.pi/2, np.pi/2, 0.0]
H = ham(k)  # Returns (n, n) numpy array

# Multiple k-points
k_list = [[0, 0, 0], [np.pi, 0, 0], [np.pi, np.pi, 0]]
H_list = ham(k_list)  # Returns list of (n, n) arrays
```

### Julia
```julia
using Firefly

# Load config
load_config!("input.cfg")

# Create Hamiltonian
ham = Hamiltonian()

# Check status
println("Uses preloaded: ", uses_preloaded(ham))

# Single k-point
k = [π/2, π/2, 0.0]
H = ham(k)  # Returns Matrix{ComplexF32}

# Multiple k-points
k_list = [[0.0, 0.0, 0.0], [π, 0.0, 0.0]]
H_list = ham(k_list)  # Returns Vector{Matrix{ComplexF32}}
```

## Config File Setup

To use preloaded Hamiltonians, set these config values:

```ini
[CONTROL]
    automatic_file_read = true
    prefix = nonexistent  # Or any name where file won't be found
    verbosity = medium

[SYSTEM]
    method = 'sc'        # or 'fcc', 'bcc', 'fermi_gas', 'emery'
    t1 = 1.0            # Nearest-neighbor hopping
    t2 = 0.2            # 2nd nearest-neighbor
    t3 = 0.1            # 3rd nearest-neighbor
    fermi_energy = 0.0  # Chemical potential
```

## Model Parameters

### Single-Band Models (SC, FCC, BCC)
- `t1`: 1st nearest-neighbor hopping
- `t2`: 2nd nearest-neighbor hopping
- `t3`: 3rd nearest-neighbor hopping
- `fermi_energy`: Chemical potential μ

### Fermi Gas
- `fermi_energy`: Chemical potential μ
- Mass is currently hardcoded to 1.0 (can be added to config)

### Emery Model (3×3)
- `t1`: t_pd (Cu-O hopping), default 1.3
- `t2`: t_pp (O-O hopping), default 0.65
- `fermi_energy`: Chemical potential μ
- ε_d and ε_p are hardcoded (can be added to config)

## Files Modified

1. `src/hamiltonian/preloaded_hamiltonians.hpp` (new)
2. `src/hamiltonian/preloaded_hamiltonians.cpp` (new)
3. `src/hamiltonian/README_PRELOADED.md` (new)
4. `src/hamiltonian/tests/test_preloaded.cpp` (new)
5. `src/hamiltonian/tests/all.cpp` (modified)
6. `src/objects/CMField/hamiltonian.hpp` (modified)
7. `src/objects/CMField/hamiltonian.cpp` (modified)
8. `src/module/exports/cpp_export.cpp` (modified)
9. `src/module/imports/cpp_imports.py` (modified)
10. `src/module/imports/cpp_imports.jl` (modified)

## Test Status

- **C++ Tests:** 14/14 passing ✓
- **Python Tests:** All tests passing ✓
- **Julia Tests:** Bindings updated and ready ✓
- **Integration Tests:** All passing ✓

## Next Steps (Optional)

Potential enhancements:
1. Add more lattice types (triangular, honeycomb, kagome)
2. Add more multi-band models (Hubbard, t-J)
3. Make Emery model parameters (ε_d, ε_p) configurable
4. Make Fermi gas mass configurable
5. Add spin-orbit coupling to models
6. Add band structure plotting utilities
