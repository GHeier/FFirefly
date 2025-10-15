#include "node.hpp"
#include "../config/load/cpp_config.hpp"
#include "../config/load/py_interface.h"
#include "band_structure.hpp"
#include "dos_tetrabz.hpp"
#include "fs.hpp"
#include <iostream>

using namespace std;

void dos_python_impl() {
    string folder = "hamiltonian";
    string filename = "dos_python";
    string function = "main";
    call_python_func(folder.c_str(), filename.c_str(), function.c_str());
}

/**
 * Wrapper function for hamiltonian category
 * Dispatches to appropriate calculation/method based on config variables
 */
extern "C" void hamiltonian_wrapper() {
    printv("Running hamiltonian_wrapper\n");
    if (calculation == "bands") {
        bands();
    }
    else if (calculation == "dos") {
        if (method == "libtetrabz") {
            dos_libtetrabz();
        }
        else if (method == "python") {
            dos_python_impl();
        }
        else
            cout << "method \"" << method << "\" not recognized for calculation dos" << endl;
    }
    else if (calculation == "fs") {
        fs();
    }
    else
        cout << "calculation \"" << calculation << "\" not recognized for category hamiltonian" << endl;
}
