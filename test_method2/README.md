# run_python_method2 and run_julia_method2 Test Suite

This directory contains tests and examples for the new shell-based execution methods.

## Quick Start

```bash
# Run standalone unit tests
./test_standalone

# Run integration tests
./integration_test.sh
```

## What's New

Added two new functions to FFirefly:
- `run_python_method2()` - Execute Python scripts via shell instead of embedded interpreter
- `run_julia_method2()` - Execute Julia scripts via shell instead of embedded interpreter

## Files

### Implementation (in main codebase)
- `src/config/load/cpp_config.hpp` - Function declarations
- `src/config/load/cpp_config.cpp` - Function implementations

### Tests
- `test_standalone.cpp` - Standalone C++ test program
- `test_standalone` - Compiled executable
- `integration_test.sh` - Shell integration tests
- `test_method2.cpp` - Alternative C++ test (not used in build)

### Test Data
- `test_config.cfg` - Sample config for unit tests
- `test_gaussian.cfg` - Config for Gaussian DOS calculation
- `test_julia_simple.jl` - Minimal Julia stdin test

### Test Scripts
- `test_category/test_calc/test_py/run.py` - Python test script
- `test_category/test_calc/test_jl/run.jl` - Julia test script

### Documentation
- `USAGE_EXAMPLE.md` - Detailed usage guide and examples

## Test Results

All tests passed successfully:

### Standalone Test
```
✓ Python method2: PASSED
✓ Julia method2:  PASSED
```

### Integration Test
```
✓ Python method2 (Gaussian DOS): PASSED
✓ Julia method2 (stdin test):     PASSED
```

## How It Works

The method2 functions use simple shell redirection:

**Python:**
```bash
python3 src/{category}/{calculation}/{method}/run.py < input.cfg
```

**Julia:**
```bash
julia src/{category}/{calculation}/{method}/run.jl < input.cfg
```

This is simpler than the original methods which use:
- `call_python_func()` with Python.h C API
- `call_julia_func()` with julia.h C API

## Advantages of Method2

1. **Simpler debugging** - Can run scripts manually
2. **No embedded interpreter complexity** - Just shell calls
3. **Process isolation** - Each run is independent
4. **Easy testing** - Scripts can be tested standalone

## Usage Example

```cpp
#include "src/config/load/cpp_config.hpp"

// In your node.cpp wrapper:
extern "C" void my_wrapper() {
    if (calculation == "test" && method == "my_python_method") {
        run_python_method2("my_python_method");
    }

    if (calculation == "test" && method == "my_julia_method") {
        run_julia_method2("my_julia_method");
    }
}
```

## See Also

- `USAGE_EXAMPLE.md` - Comprehensive usage guide
- Main FFirefly documentation: `CLAUDE.md`
