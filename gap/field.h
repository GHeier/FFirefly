#pragma once
#ifndef field_class_H
#define field_class_H

#include <string>
#include "vec.h"

using namespace std;

class ScalarField {
    public:
        vector<Vec> points;
        vector<float> values;
        float xmin, xmax, ymin, ymax, zmin, zmax, wmin, wmax;
        int nx, ny, nz, nw;
        int dimension;

        ScalarField();
        ScalarField(vector<Vec> points, vector<float> values);
        ScalarField(string filename, int dimension=3);

        void get_values_for_interpolation();
        float operator() (Vec point);
};

class VectorField {
    public:
        vector<Vec> points;
        vector<Vec> values;
        float xmin, xmax, ymin, ymax, zmin, zmax, wmin, wmax;
        int nx, ny, nz, nw;
        int dimension;

        VectorField();
        VectorField(vector<Vec> points, vector<Vec> values);
        VectorField(string filename, int dimension=3);

        void get_values_for_interpolation();
        Vec operator() (Vec point);
};

#endif
