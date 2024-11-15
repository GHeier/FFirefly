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
        float xmin, xmax, ymin, ymax, zmin, zmax, wmin, wmax;
        int nx, ny, nz, nw;
        int dimension;
        bool is_complex;

        ScalarField();
        ScalarField(vector<Vec> points, vector<float> values);
        ScalarField(vector<Vec> points, vector<complex<float>> values);
        ScalarField(string filename, int dimension=3, bool is_complex=false);

        void get_values_for_interpolation();
        float operator() (Vec point);
};

class VectorField {
    public:
        vector<Vec> points;
        vector<Vec> values;
        vector<Vec> ivalues;
        float xmin, xmax, ymin, ymax, zmin, zmax, wmin, wmax;
        int nx, ny, nz, nw;
        int dimension;
        bool is_complex;

        VectorField();
        VectorField(vector<Vec> points, vector<Vec> values);
        VectorField(vector<Vec> points, vector<complex<float>> values);
        VectorField(string filename, int dimension=3, bool is_complex=false);

        void get_values_for_interpolation();
        Vec operator() (Vec point);
};

#endif
