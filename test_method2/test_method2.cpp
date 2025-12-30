#include "cpp_config.hpp"
#include <iostream>
#include <fstream>
#include <cstdio>

using namespace std;

// Helper to create a test config file
void create_test_config(const string& filename,
                       const string& cat,
                       const string& calc,
                       const string& meth) {
    ofstream f(filename);
    f << "[CONTROL]\n";
    f << "    category = '" << cat << "'\n";
    f << "    calculation = '" << calc << "'\n";
    f << "    method = '" << meth << "'\n";
    f << "    prefix = 'test'\n";
    f << "    verbosity = 'high'\n";
    f << "\n";
    f << "[SYSTEM]\n";
    f << "    nbnd = 1\n";
    f.close();
}

int main() {
    cout << "=== Testing run_python_method2 and run_julia_method2 ===" << endl;

    // Test 1: Python method2
    cout << "\n--- Test 1: Python method2 ---" << endl;
    create_test_config("test_py_method2.cfg", "hamiltonian", "DOS", "gaussian");

    // Load config
    read_c_config_wrapper("test_py_method2.cfg");

    cout << "Category: " << category << endl;
    cout << "Calculation: " << calculation << endl;
    cout << "Method: gaussian" << endl;

    // Copy config to input.cfg (required by method2)
    system("cp test_py_method2.cfg input.cfg");

    int result_py = run_python_method2("gaussian");
    cout << "Python method2 exit code: " << result_py << endl;

    if (result_py == 0) {
        cout << "✓ Python method2 test PASSED" << endl;
    } else {
        cout << "✗ Python method2 test FAILED" << endl;
    }

    // Test 2: Julia method2
    cout << "\n--- Test 2: Julia method2 ---" << endl;
    create_test_config("test_jl_method2.cfg", "superconductor", "eliashberg", "power_iteration");

    // Load config
    read_c_config_wrapper("test_jl_method2.cfg");

    cout << "Category: " << category << endl;
    cout << "Calculation: " << calculation << endl;
    cout << "Method: power_iteration" << endl;

    // Copy config to input.cfg (required by method2)
    system("cp test_jl_method2.cfg input.cfg");

    int result_jl = run_julia_method2("power_iteration");
    cout << "Julia method2 exit code: " << result_jl << endl;

    if (result_jl == 0) {
        cout << "✓ Julia method2 test PASSED" << endl;
    } else {
        cout << "✗ Julia method2 test FAILED (may need proper input files)" << endl;
    }

    // Cleanup
    remove("test_py_method2.cfg");
    remove("test_jl_method2.cfg");
    remove("input.cfg");

    cout << "\n=== Test Summary ===" << endl;
    if (result_py == 0) {
        cout << "Python method2: PASSED" << endl;
    } else {
        cout << "Python method2: FAILED" << endl;
    }

    if (result_jl == 0) {
        cout << "Julia method2: PASSED" << endl;
    } else {
        cout << "Julia method2: FAILED (this is expected if input files are missing)" << endl;
    }

    return 0;
}
