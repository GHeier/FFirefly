# Category Generator Usage

The `write.py` script automatically generates the directory structure and boilerplate files for FFirefly categories.

## Fixed Issues

1. **Missing import**: Added `import shutil` at the top
2. **Path handling**: Converted all string concatenation to proper `Path` objects
3. **Template directory**: Fixed relative paths to use `Path(__file__).parent` for template location
4. **Function signatures**: Added `template_dir` parameter to `make_run_file` and `make_test_file`
5. **Directory management**: Fixed the nested loop that was incorrectly mutating `dir` variable
6. **Node file generation**: Fixed to use correct parameters (removed unused `methods` parameter, renamed `dir` to `category_dir`)
7. **Test file location**: Fixed to create test files in the `tests/` subdirectory
8. **Missing template**: Created `hpp_test.txt` template file
9. **Return value**: Made `make_folders` return useful records of what was created

## Usage

### As a Python Module

```python
from pathlib import Path
from write import make_folders

categories = {
    "hamiltonian": {
        "generate": {
            "hr_from_hk": "python"
        },
        "DOS": {
            "tetrahedra": "c++"
        },
        "FS": {
            "tetrahedra": "c++"
        }
    },
    "superconductor": {
        "eliashberg": {
            "power_iteration": "julia"
        }
    }
}

# Create folders in default location (../../ relative to script)
records = make_folders(categories)

# Or specify custom base directory
records = make_folders(categories, base_dir="/custom/path")
```

### As a Standalone Script

```bash
cd /home/g/Research/FFirefly/src/config/categories
python write.py
```

This runs the example at the bottom of the file.

## What Gets Created

For each `category/calculation/method` combination:

```
src/
  category/
    node.cpp              # Auto-generated category wrapper
    node.hpp              # Header for category wrapper
    shared/               # Shared code for category
    calculation/
      method/
        README.md         # Documentation template
        run.{ext}         # From template (cpp/py/jl)
        run.hpp           # For C++ only
        tests/
          test.{ext}      # From template
          test.hpp        # For C++ only
```

## Template Files

The script copies from these templates in `src/config/categories/`:
- `README_template.txt` → `README.md` (with {CATEGORY}/{CALCULATION}/{METHOD} substitutions)
- `cpp_run.txt` → `run.cpp`
- `hpp_run.txt` → `run.hpp`
- `cpp_test.txt` → `test.cpp`
- `hpp_test.txt` → `test.hpp`
- `python_run.txt` → `run.py`
- `python_test.txt` → `test.py`
- `julia_run.txt` → `run.jl`
- `julia_test.txt` → `test.jl`

### README Template Structure

Each method gets a `README.md` with these sections:
- **Overview** - Purpose in FFirefly framework
- **Quick Description** - One-line method summary
- **Dependencies** - Required and optional dependencies
- **Install Instructions** - Setup steps if needed
- **Results Saved** - Output file descriptions
- **Testing** - How to run tests
- **Calculation Details** - Algorithm, parameters, implementation notes
- **References** - Citations and links

## Generated Node Files

The `node.cpp` file is auto-generated with proper if/else chains:

```cpp
extern "C" void hamiltonian_wrapper() {
    printv("Running hamiltonian_wrapper\n");
    if (calculation == "generate" && method == "hr_from_hk") run_python_method("hr_from_hk");
    else if (calculation == "DOS" && method == "tetrahedra") run_cpp_method("tetrahedra");
    else if (calculation == "FS" && method == "tetrahedra") run_cpp_method("tetrahedra");
    else {
        printf("In hamiltonian category, calculation `%s` with method `%s` not recognized\n",
               calculation.c_str(), method.c_str());
    }
}
```

## Safe Overwrites

The script only creates files if they don't exist:
- Existing `README.md` files are preserved
- Existing `run.*` files are preserved
- Existing `test.*` files are preserved
- Directories are created with `exist_ok=True`
- Node files are always regenerated (to keep dispatch logic in sync)

This allows you to:
- Re-run the script safely to add new methods
- Regenerate node files after adding/removing methods
- Preserve your custom implementations and documentation
