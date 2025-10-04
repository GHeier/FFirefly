#include "hamiltonian.hpp"
#include "../../config/load/cpp_config.hpp"
#include <complex>
#include <fstream>
#include <iostream>

using namespace std;

Hamiltonian::Hamiltonian() {
    string fname = indir + prefix + "_hamiltonian.h5";
    file_found = false;

    // Try to load from file if automatic_file_read is enabled
    if (automatic_file_read) {
        ifstream f(fname.c_str());
        if (f.good()) {
            field = Field_CM(fname.c_str());
            file_found = true;
            if (verbosity == "high")
                cout << "Loaded Hamiltonian from " << fname << endl;
        } else {
            if (verbosity == "high")
                cout << "Hamiltonian file not found: " << fname << endl;
        }
    }
}

vector<vector<complex<float>>> Hamiltonian::operator()(Vec k, float w) {
    if (!file_found) {
        cerr << "Error: Hamiltonian not loaded from file" << endl;
        // Return empty matrix
        return vector<vector<complex<float>>>();
    }
    return field(k, w);
}
