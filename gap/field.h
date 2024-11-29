#pragma once
#ifndef field_class_H
#define field_class_H

#include <string>
#include <complex>
#include "vec.h"

using namespace std;

class ScalarField {
    public:
        vector<Vec> points;
        vector<float> values;
        vector<float> ivalues;
        vector<Vec> domain;
        vector<Vec> inv_domain;
        float xmin, xmax, ymin, ymax, zmin, zmax, wmin, wmax;
        int nx, ny, nz, nw;
        int dimension;
        bool is_complex;

        ScalarField();
        ScalarField(vector<Vec> points, vector<float> values, int dimension=3, bool is_complex=false);
        ScalarField(vector<Vec> points, vector<complex<float>> values, int dimension=3, bool is_complex=true);
        ScalarField(string filename, int dimension=3, bool is_complex=false);

        void get_values_for_interpolation();
        float operator() (Vec point);
};

class VectorField {
    public:
        vector<Vec> points;
        vector<Vec> values;
        vector<Vec> ivalues;
        vector<Vec> domain;
        vector<Vec> inv_domain;
        float xmin, xmax, ymin, ymax, zmin, zmax, wmin, wmax;
        int nx, ny, nz, nw;
        int dimension;
        bool is_complex;

        VectorField();
        VectorField(vector<Vec> points, vector<Vec> values, int dimension=3, bool is_complex=false);
        VectorField(vector<Vec> points, vector<complex<float>> values, int dimension=3, bool is_complex=true);
        VectorField(string filename, int dimension=3, bool is_complex=false);

        void get_values_for_interpolation();
        Vec operator() (Vec point);
};

bool domain_vec_found(float a[4], float dx, float dy, float dz, float dw);
vector<Vec> invertMatrix(vector<Vec>& matrix, int n);

#endif
