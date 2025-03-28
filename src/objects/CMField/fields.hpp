/*
 * This file contains the definition of the ScalarField and VectorField classes.
 * These classes are used to represent scalar and vector fields, respectively.
 * They are implemented over a mesh of points and values, and can be used to interpolate
 * The mesh is interpolated automatically when the field is created. Ideal for BZ calculations
 *
 * Author: Griffin Heier
 */
#pragma once

#include <complex>

#include "../vec.hpp"
#include "cmfield.hpp"

using namespace std;
class Field_C {
    public:
        CMF cmf;

        Field_C();
        Field_C(CMF cmf);
        Field_C(string filename);

        complex<float> operator() (Vec point, float w = 0);
        complex<float> operator() (float w);
};

class Field_R {
    public:
        CMF cmf;

        Field_R();
        Field_R(CMF cmf);
        Field_R(string filename);

        float operator() (Vec point, float w = 0);
        float operator() (double w);
};

extern "C" float test_func(float w);
