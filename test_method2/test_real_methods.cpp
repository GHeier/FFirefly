#include "../src/config/load/cpp_config.hpp"
#include <iostream>

using namespace std;

extern "C" void load_cpp_config();

int main() {
    cout << "=== Testing real run_python_method2 and run_julia_method2 ===" << endl;

    // Test 1: Python method2 with gaussian DOS
    cout << "\n--- Test 1: Python method2 (gaussian DOS) ---" << endl;

    // Load the config
    read_c_config_wrapper("test_gaussian.cfg");

    cout << "Config loaded:" << endl;
    cout << "  Category: " << category << endl;
    cout << "  Calculation: " << calculation << endl;
    cout << "  Method: gaussian" << endl;

    // Copy to input.cfg for method2
    system("cp test_gaussian.cfg input.cfg");

    // Test run_python_method2
    int result = run_python_method2("gaussian");
    cout << "\nPython method2 exit code: " << result << endl;

    if (result == 0) {
        cout << "✓ Python method2 (gaussian) PASSED" << endl;
    } else {
        cout << "✗ Python method2 (gaussian) FAILED with exit code " << result << endl;
    }

    // Cleanup
    system("rm -f input.cfg");

    cout << "\n=== Test Complete ===" << endl;

    return result;
}
