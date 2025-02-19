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

#include "vec.hpp"
#include "cmf.hpp"

using namespace std;
class CMF_CS {
    public:
        CMF base;

        CMF_CS();
        CMF_CS(CMF base);
        CMF_CS(string filename);

        complex<float> operator() (Vec point, float w = 0);
        complex<float> operator() (float w);
};

class CMF_RS {
    public:
        CMF base;

        CMF_RS();
        CMF_RS(CMF base);
        CMF_RS(string filename);

        float operator() (Vec point, float w = 0);
        float operator() (float w);
};

