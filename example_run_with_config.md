# run_with_config Function

A new C++ function has been added to `src/config/load/cpp_config.cpp` to run executables with config file input via stdin.

## Function Signature

```cpp
int run_with_config(const std::string& executable, const std::string& config_file);
```

## Description

Runs an executable with a configuration file piped to stdin, equivalent to the shell command:
```bash
executable < config_file
```

## Parameters

- `executable` - Path to the executable file (e.g., `"./program.exe"`, `"build/bin/fly.x"`)
- `config_file` - Path to the configuration file (e.g., `"input.cfg"`)

## Return Value

- Returns the exit code of the executable (0 typically means success)
- Returns -1 if:
  - The executable file doesn't exist
  - The config file doesn't exist
  - The command fails to execute
  - The process terminates abnormally

## Features

- **File validation**: Checks that both executable and config file exist before running
- **Verbosity support**: Shows the command being run when `verbosity = "high"`
- **Proper exit code handling**: Uses `WEXITSTATUS` to extract the actual exit code
- **Error reporting**: Provides clear error messages for common failure cases

## Usage Example

```cpp
#include "src/config/load/cpp_config.hpp"

int main() {
    // Run fly.x with input.cfg
    int result = run_with_config("build/bin/fly.x", "input.cfg");

    if (result == 0) {
        std::cout << "Execution successful!\n";
    } else if (result == -1) {
        std::cerr << "Error: Failed to run executable\n";
    } else {
        std::cerr << "Executable exited with code: " << result << "\n";
    }

    return result;
}
```

## Implementation Details

The function:
1. Uses `std::filesystem::exists()` to validate paths
2. Constructs the shell command with input redirection
3. Executes via `std::system()`
4. Extracts the proper exit code using `WIFEXITED()` and `WEXITSTATUS()` macros

## Headers Required

```cpp
#include "src/config/load/cpp_config.hpp"
```

The implementation uses:
- `<cstdlib>` for `std::system()`
- `<sys/wait.h>` for `WIFEXITED()` and `WEXITSTATUS()`
- `<filesystem>` for file existence checks
