# Using run_python_method2 and run_julia_method2

## Overview

The new `run_python_method2` and `run_julia_method2` functions provide an alternative way to execute Python and Julia scripts using **pure shell interaction** instead of embedded interpreters (Python.h / julia.h).

## Key Differences

| Feature | Original Methods | Method2 (New) |
|---------|-----------------|---------------|
| **Execution** | Uses embedded Python/Julia via C API | Uses shell system() calls |
| **Dependencies** | Requires Python.h and julia.h headers | No special headers needed |
| **Process** | In-process execution | Separate process via shell |
| **Config Input** | Via run_with_config() helper | Direct stdin redirection (`< input.cfg`) |
| **Startup Cost** | Reuses interpreter | New interpreter per call |

## Implementation

Both methods are implemented in:
- **Header**: `src/config/load/cpp_config.hpp`
- **Source**: `src/config/load/cpp_config.cpp`

```cpp
// Method signatures
int run_python_method2(const std::string& method_name);
int run_julia_method2(const std::string& method_name);
```

## How They Work

### run_python_method2

```cpp
int run_python_method2(const string& method_name) {
    // 1. Locate the script: src/{category}/{calculation}/{method_name}/run.py
    // 2. Execute: python3 script.py < input.cfg
    // 3. Return exit code
}
```

### run_julia_method2

```cpp
int run_julia_method2(const string& method_name) {
    // 1. Locate the script: src/{category}/{calculation}/{method_name}/run.jl
    // 2. Execute: julia script.jl < input.cfg
    // 3. Return exit code
}
```

## Usage in node.cpp Files

You can use method2 functions exactly like the original methods:

```cpp
// Example: src/hamiltonian/node.cpp
extern "C" void hamiltonian_wrapper() {
    printv("Running hamiltonian_wrapper\n");

    // Original method (embedded Python)
    if (calculation == "DOS" && method == "gaussian")
        run_python_method("gaussian");

    // Alternative: method2 (shell-based)
    if (calculation == "DOS" && method == "gaussian_v2")
        run_python_method2("gaussian");

    // For Julia
    if (calculation == "test" && method == "julia_test")
        run_julia_method2("julia_test");
}
```

## When to Use Method2

### Advantages
- ✓ Simpler debugging (can run scripts manually)
- ✓ No dependency on Python.h/julia.h internals
- ✓ Clean process isolation
- ✓ Easy to test scripts independently

### Disadvantages
- ✗ Higher startup cost (new interpreter each time)
- ✗ Cannot share state between calls
- ✗ Slightly slower for repeated calls

### Recommended Use Cases
1. **Development/Testing**: Quick iteration on Python/Julia scripts
2. **One-shot calculations**: Methods that run once per session
3. **Debugging**: When embedded interpreter causes issues
4. **Isolation**: When you want clean process separation

## Testing

### Unit Tests
Run the standalone test:
```bash
cd test_method2
./test_standalone
```

### Integration Tests
```bash
cd test_method2
./integration_test.sh
```

### Manual Testing
```bash
# Test Python method2 directly
cd /home/g/Research/FFirefly
python3 src/hamiltonian/DOS/gaussian/run.py < test_method2/test_gaussian.cfg

# Test Julia method2 directly
julia test_method2/test_julia_simple.jl < test_method2/test_gaussian.cfg
```

## Example Scripts

See the test directory for examples:
- `test_method2/test_category/test_calc/test_py/run.py` - Simple Python test
- `test_method2/test_category/test_calc/test_jl/run.jl` - Simple Julia test
- `test_method2/test_julia_simple.jl` - Minimal Julia stdin test

## Notes

1. **Config File**: Both methods expect `input.cfg` in the current directory (created by run_with_config internally for original methods, but method2 uses direct stdin redirection)

2. **Exit Codes**: Both methods return the script's exit code (0 = success, non-zero = failure)

3. **Verbosity**: If `verbosity = 'high'` in config, the command being executed will be printed

4. **Error Handling**: Methods return -1 if:
   - Cannot find src directory
   - Script execution fails
   - Process doesn't terminate normally

5. **Script Location**: Scripts must follow the standard path convention:
   ```
   src/{category}/{calculation}/{method_name}/run.{py,jl}
   ```
