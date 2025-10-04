#pragma once
#include "fields.hpp"

#include <complex>
#include <vector>

using namespace std;

class Hamiltonian {
  public:
    Field_CM field;
    bool file_found;

    Hamiltonian();
    // Returns matrix H_ab(k) where result[a][b] is the Hamiltonian element
    vector<vector<complex<float>>> operator()(Vec k, float w = 0);
};
