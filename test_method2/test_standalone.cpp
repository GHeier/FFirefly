#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <sys/wait.h>
#include <unistd.h>
#include <linux/limits.h>

using namespace std;

// Simplified versions of the method2 functions for standalone testing
string get_test_loc() {
    char path[PATH_MAX];
    ssize_t len = readlink("/proc/self/exe", path, sizeof(path) - 1);
    if (len == -1) {
        return "./";
    }
    path[len] = '\0';

    string exe_path(path);
    size_t last_slash = exe_path.find_last_of('/');
    if (last_slash != string::npos) {
        return exe_path.substr(0, last_slash + 1);
    }
    return "./";
}

int run_python_method2_test(const string& category, const string& calculation, const string& method_name) {
    string loc = get_test_loc();
    string script_path = loc + category + "/" + calculation + "/" + method_name + "/run.py";

    string command = "python3 " + script_path + " < input.cfg";

    cout << "Running (method2): " << command << endl;

    int result = system(command.c_str());

    if (result == -1) {
        cerr << "Error: Failed to execute Python script (method2).\n";
        return -1;
    }

    if (WIFEXITED(result)) {
        return WEXITSTATUS(result);
    } else {
        cerr << "Error: Python process did not terminate normally (method2).\n";
        return -1;
    }
}

int run_julia_method2_test(const string& category, const string& calculation, const string& method_name) {
    string loc = get_test_loc();
    string script_path = loc + category + "/" + calculation + "/" + method_name + "/run.jl";

    string command = "julia " + script_path + " < input.cfg";

    cout << "Running (method2): " << command << endl;

    int result = system(command.c_str());

    if (result == -1) {
        cerr << "Error: Failed to execute Julia script (method2).\n";
        return -1;
    }

    if (WIFEXITED(result)) {
        return WEXITSTATUS(result);
    } else {
        cerr << "Error: Julia process did not terminate normally (method2).\n";
        return -1;
    }
}

int main() {
    cout << "========================================" << endl;
    cout << "Testing run_python_method2 and run_julia_method2" << endl;
    cout << "========================================" << endl;

    // Create input.cfg
    system("cp test_config.cfg input.cfg");

    // Test 1: Python method2
    cout << "\n--- Test 1: run_python_method2 ---" << endl;
    int result_py = run_python_method2_test("test_category", "test_calc", "test_py");
    cout << "Exit code: " << result_py << endl;

    if (result_py == 0) {
        cout << "✓ Python method2 test PASSED" << endl;
    } else {
        cout << "✗ Python method2 test FAILED" << endl;
    }

    // Test 2: Julia method2
    cout << "\n--- Test 2: run_julia_method2 ---" << endl;
    int result_jl = run_julia_method2_test("test_category", "test_calc", "test_jl");
    cout << "Exit code: " << result_jl << endl;

    if (result_jl == 0) {
        cout << "✓ Julia method2 test PASSED" << endl;
    } else {
        cout << "✗ Julia method2 test FAILED" << endl;
    }

    // Summary
    cout << "\n========================================" << endl;
    cout << "Summary:" << endl;
    cout << "  Python method2: " << (result_py == 0 ? "PASSED ✓" : "FAILED ✗") << endl;
    cout << "  Julia method2:  " << (result_jl == 0 ? "PASSED ✓" : "FAILED ✗") << endl;
    cout << "========================================" << endl;

    // Cleanup
    remove("input.cfg");

    return (result_py == 0 && result_jl == 0) ? 0 : 1;
}
