# Category Code Generation System

This directory contains automatic code generators for FFirefly's category system.

## Overview

The category system automatically generates:
- **node.cpp/node.hpp** files for each category
- **main.c** integration code (includes, functions, dispatchers)

This ensures consistency across the codebase and reduces boilerplate when adding new categories or calculations.

## Files

### `categories.py`
Master definition file containing:
- `CATEGORIES`: Dict mapping categories to their calculations and implementation files
- `CATEGORIES_WITH_TESTS`: List of categories that have test suites

**File extension determines calling convention:**
- `.cpp` → Direct C++ function call
- `.py` → `call_python_func()` wrapper
- `.jl` → `call_julia_func()` wrapper with module name

### `write_nodes.py`
Generates `node.cpp` and `node.hpp` files for each category:
- Creates `extern "C"` wrapper functions
- Includes appropriate interfaces (Python/Julia/C++)
- Dispatches based on `calculation` config variable
- Generates helper functions for Python/Julia calculations

### `write_main.py`
Generates main.c integration code:
- Include statements for all category nodes
- Function definitions that call category wrappers
- Dispatcher (if-else chain) for category routing
- Test function with all test calls

## Usage

### Running the Generators

```bash
cd /home/g/Research/FFirefly/src
python config/categories.py
```

This creates files in `config/categories/archive/`:
- `main_sections.txt` - Code to integrate into main.c
- `<category>/node.cpp` - Category wrapper implementations
- `<category>/node.hpp` - Category wrapper declarations

### Adding a New Category

1. Edit `config/categories.py`:
```python
CATEGORIES = {
    "my_category": {
        "calculation1": "file1.cpp",
        "calculation2": "file2.py",
        "calculation3": "file3.jl",
    },
}
```

2. Run the generator:
```bash
python config/categories.py
```

3. Copy generated files to your category directory:
```bash
cp config/categories/archive/my_category/node.* src/my_category/
```

4. Update `main.c` with sections from `archive/main_sections.txt`

### Adding a New Calculation

1. Add entry to existing category in `categories.py`:
```python
"hamiltonian": {
    "fs": "fs.cpp",
    "dos": "dos.cpp",
    "my_new_calc": "my_new_calc.py",  # Add this
},
```

2. Regenerate nodes:
```bash
python config/categories.py
```

3. Replace old node files with new ones

## Generated Code Structure

### node.hpp
```cpp
#pragma once

#ifdef __cplusplus
extern "C" {
#endif

void category_wrapper();

#ifdef __cplusplus
}
#endif
```

### node.cpp
```cpp
#include "node.hpp"
#include "../config/load/cpp_config.hpp"
// Includes for Python/Julia if needed
#include "calculation_files.hpp"

// Helper functions for Python/Julia calls
void calc_impl() {
    call_python_func(...);
}

// Main dispatcher
extern "C" void category_wrapper() {
    if (calculation == "calc1")
        calc1();  // Direct C++ call
    else if (calculation == "calc2")
        calc2_impl();  // Python/Julia wrapper
    // ...
}
```

### main.c sections
```c
// Includes
#include "category/node.hpp"

// Function definition
void category() {
    printf("Starting Category Calculation\\n\\n");
    category_wrapper();
}

// Dispatcher
if (!strcmp(category, "category_name"))
    category();
```

## Examples

### C++ Calculation
```python
"hamiltonian": {
    "fs": "fs.cpp",  # Calls fs() directly
}
```

### Python Calculation
```python
"analysis": {
    "plot": "plot_data.py",  # Calls call_python_func("analysis", "plot_data", "main")
}
```

### Julia Calculation
```python
"many_body": {
    "dmft": "dmft_loop.jl",  # Calls call_julia_func("many_body/", "dmft_loop", "DmftLoop", "main")
}
```

## Notes

- Module names for Julia are auto-generated from filename in TitleCase
  - `my_function.jl` → `MyFunction` module
- All Python/Julia calculations must have a `main()` function
- C++ calculations can be any function name matching the filename
- Category names with slashes (e.g., "config/load") work correctly
