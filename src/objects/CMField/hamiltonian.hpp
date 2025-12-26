#pragma once
#include "fields.hpp"

#include <complex>
#include <vector>
#include <string>

using namespace std;

class Hamiltonian {
  public:
    Field_CM field;
    bool file_found;

    Hamiltonian();
    vector<vector<complex<float>>> operator()(Vec k);
    vector<vector<vector<complex<float>>>> operator()(vector<Vec> kpoints);

    vector<float> get_bands(Vec k);
    vector<vector<float>> get_bands(vector<Vec> kpoints);

    vector<eigvec> get_wavefunctions(Vec k);
    vector<vector<eigvec>> get_wavefunctions(vector<Vec> kpoints);
};

