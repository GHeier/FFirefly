#include "node.hpp"
#include "../config/load/py_interface.h"
#include "../config/load/cpp_config.hpp"
#include <iostream>
#include <string>

using namespace std;

extern "C" void hamiltonian_wrapper() {
    printv("Running hamiltonian_wrapper\n");

    if (calculation == "generate") {
        string folder = "hamiltonian/";
        string filename = "Hk";
        string function = "generate_hamiltonian";
        call_python_func(folder.c_str(), filename.c_str(), function.c_str());
    }
    else {
        cout << "calculation " << calculation << " not recognized for hamiltonian category" << endl;
        cout << "Available calculations: generate" << endl;
    }
}
