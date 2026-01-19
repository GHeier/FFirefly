#include "hamiltonian.hpp"
#include "src/hamiltonian/models/preloaded_hamiltonians.hpp"
#include "src/config/load/cpp_config.hpp"
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

vector<Vec> Hamiltonian::get_fermi_velocity(Vec k) {
    // Compute Fermi velocity v_n(k) = ∇_k E_n(k) using finite differences
    // Returns a vector of Vec, one for each band

    float dk = 0.001;  // Small step for numerical derivative
    int dim = k.dimension;

    // Get bands at current k
    vector<float> E0 = get_bands(k);
    int nbands = E0.size();

    vector<Vec> velocities(nbands);

    // Compute gradient for each band
    for (int n = 0; n < nbands; n++) {
        Vec v;
        v.dimension = dim;

        // dx derivative
        Vec kx_plus = k;
        kx_plus.x += dk;
        vector<float> Ex_plus = get_bands(kx_plus);
        v.x = (Ex_plus[n] - E0[n]) / dk;

        // dy derivative
        Vec ky_plus = k;
        ky_plus.y += dk;
        vector<float> Ey_plus = get_bands(ky_plus);
        v.y = (Ey_plus[n] - E0[n]) / dk;

        // dz derivative (only for 3D)
        if (dim == 3) {
            Vec kz_plus = k;
            kz_plus.z += dk;
            vector<float> Ez_plus = get_bands(kz_plus);
            v.z = (Ez_plus[n] - E0[n]) / dk;
        } else {
            v.z = 0.0;
        }

        velocities[n] = v;
    }

    return velocities;
}

vector<vector<Vec>> Hamiltonian::get_fermi_velocity(vector<Vec> kpoints) {
    vector<vector<Vec>> velocities;
    for (auto k : kpoints) {
        velocities.push_back(get_fermi_velocity(k));
    }
    return velocities;
}
