#include "hamiltonian.hpp"
#include "../../hamiltonian/preloaded_hamiltonians.hpp"
#include "../../config/load/cpp_config.hpp"
#include "fields.hpp"
#include <complex>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cctype>

using namespace std;

Hamiltonian::Hamiltonian() {
    string fname = outdir + prefix + "_hamiltonian.h5";
    file_found = false;

    ifstream file(fname.c_str());
    if (file.good()) {
        file_found = true;
        printv("Loading Hamiltonian from file: %s\n\n", fname.c_str());
        field = Field_CM(fname.c_str());
    }
    else {
        printv("Hamiltonian file not found: %s\n", fname.c_str());
        printv("Using preloaded hamiltonian model: %s\n\n", hamiltonian.c_str());
    }
}

vector<vector<complex<float>>> Hamiltonian::operator()(Vec k) {
    if (file_found)
        return field(k, 0);
    else 
        return H(k);
}

vector<vector<vector<complex<float>>>> Hamiltonian::operator()(vector<Vec> kpoints) {
    vector<vector<vector<complex<float>>>> Hkpoints;
    for (auto k : kpoints) {
        Hkpoints.push_back(operator()(k));
    }
    return Hkpoints;
}

vector<float> Hamiltonian::get_bands(Vec k) {
    vector<vector<complex<float>>> Hk = operator()(k);
    return diag(Hk);
}

vector<vector<float>> Hamiltonian::get_bands(vector<Vec> kpoints) {
    vector<vector<float>> bands;
    for (auto k : kpoints) {
        bands.push_back(get_bands(k));
    }
    return bands;
}

vector<eigvec> Hamiltonian::get_wavefunctions(Vec k) {
    vector<vector<complex<float>>> Hk = operator()(k);
    return fulldiag(Hk);
}

vector<vector<eigvec>> Hamiltonian::get_wavefunctions(vector<Vec> kpoints) {
    vector<vector<eigvec>> wavefunctions;
    for (auto k : kpoints) {
        wavefunctions.push_back(get_wavefunctions(k));
    }
    return wavefunctions;
}
