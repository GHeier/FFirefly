# Typed Language Interface Functions

This document summarizes the typed interface functions for calling Python and Julia from C/C++.

## Overview

The codebase now supports calling Python and Julia functions with typed return values:
- **bool** - Boolean values (true/false)
- **int** - 32-bit integers
- **float** - 32-bit floating point
- **double** - 64-bit floating point
- **string** - C-strings (char*)

## Python Interface Functions

### Location
- **Header**: `src/config/load/py_interface.h`
- **Implementation**: `src/config/load/py_interface.c`

### Functions

```c
bool call_python_func_bool(const char *folder, const char *filename, const char *function);
int call_python_func_int(const char *folder, const char *filename, const char *function);
float call_python_func_float(const char *folder, const char *filename, const char *function);
double call_python_func_double(const char *folder, const char *filename, const char *function);
const char* call_python_func_string(const char *folder, const char *filename, const char *function);
```

### Usage Example

```c
#include "config/load/py_interface.h"

start_python();

// Call Python function that returns int
int result = call_python_func_int("config", "test_typed_interface", "return_int_42");
// result = 42

// Call Python function that returns float
float pi = call_python_func_float("config", "test_typed_interface", "return_float_pi");
// pi = 3.14159

// Call Python function that returns string
const char* msg = call_python_func_string("config", "test_typed_interface", "return_string_hello");
// msg = "Hello, World!"
free((void*)msg);  // Must free the returned string!

end_python();
```

## Julia Interface Functions

### Location
- **Header**: `src/config/load/jl_interface.h`
- **Implementation**: `src/config/load/jl_interface.c`

### Functions

```c
bool call_julia_func_bool(const char *folder, const char *filename, const char *module, const char *function);
int call_julia_func_int(const char *folder, const char *filename, const char *module, const char *function);
float call_julia_func_float(const char *folder, const char *filename, const char *module, const char *function);
double call_julia_func_double(const char *folder, const char *filename, const char *module, const char *function);
const char* call_julia_func_string(const char *folder, const char *filename, const char *module, const char *function);
```

### Usage Example

```c
#include "config/load/jl_interface.h"

// Call Julia function that returns int (note the extra module parameter)
int result = call_julia_func_int("config/", "test_typed_interface", "TestTypedInterface", "return_int_42");
// result = 42

// Call Julia function that returns float
float pi = call_julia_func_float("config/", "test_typed_interface", "TestTypedInterface", "return_float_pi");
// pi = 3.14159

// Call Julia function that returns string
const char* msg = call_julia_func_string("config/", "test_typed_interface", "TestTypedInterface", "return_string_hello");
// msg = "Hello, World!"
free((void*)msg);  // Must free the returned string!
```

## Test Files

### Python Tests
- **Test Module**: `src/config/test_typed_interface.py`
- **C++ Tests**: `src/config/load/tests/py_interface_typed_tests.cpp`

### Julia Tests
- **Test Module**: `src/config/test_typed_interface.jl`
- **C++ Tests**: `src/config/load/tests/jl_interface_typed_tests.cpp`

## Running Tests

Tests are automatically run when executing `fly.x` with no arguments:

```bash
cd /path/to/anywhere
/home/g/Research/FFirefly/build/bin/fly.x
```

Expected output:
```
Running Config Load tests
✓ All 6 Python Interface tests passed!
✓ All 14 Python Typed Interface tests passed!
✓ All 16 Julia Typed Interface tests passed!
```

## Test Coverage

### Python Interface Tests (14 tests)
- ✓ Integer return: 42, -100, 0, computed sum (60)
- ✓ Float return: pi, negative, zero, computed product (10.0)
- ✓ Double return: large value, small value
- ✓ String return: "Hello, World!", empty string, special characters, computed message

### Julia Interface Tests (16 tests)
- ✓ Integer return: 42, -100, 0, computed sum (60)
- ✓ Float return: pi, negative, zero, computed product (10.0)
- ✓ Double return: large value, small value
- ✓ String return: "Hello, World!", empty string, special characters, computed message
- ✓ Boolean return: true, false

## Important Notes

1. **String Memory Management**: Returned strings from `call_python_func_string()` and `call_julia_func_string()` are allocated with `strdup()` and **must be freed** by the caller using `free()`.

2. **Python Initialization**: Must call `start_python()` before any Python interface functions and `end_python()` when done.

3. **Julia Initialization**: Each Julia function call initializes and cleans up Julia automatically via `jl_init()` and `jl_atexit_hook()`.

4. **Path Parameters**:
   - Python: `folder` is relative to `src/` (e.g., "config")
   - Julia: `folder` usually includes trailing slash (e.g., "config/")

5. **Error Handling**: All functions return sensible defaults on error:
   - `bool`: false
   - `int`: 0
   - `float`/`double`: 0.0
   - `string`: empty string ""

6. **Type Conversion**: Python functions handle automatic type conversion (e.g., int → float), but Julia functions are more strict about types.
