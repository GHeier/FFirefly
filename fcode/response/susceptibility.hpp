#pragma once

#include <vector>
#include <complex>
#include <functional>
#include "../objects/vec.hpp"
#include "../objects/field.hpp"

using namespace std;

class Susceptibility {
    public:
        Field chi;

        Susceptibility();
        Susceptibility(vector<Vec> points, vector<complex<float>> values, int dimension=3, bool is_complex=true);
        Susceptibility(string filename, int dimension=3, bool is_complex=false);

        float operator() (Vec point, float w);
};


float fermi_dirac(float E, float T);

float ratio(Vec k, Vec q, float w, float T);

Vec to_IBZ(const Vec k);

