#pragma once

#include "vec.hpp"

#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace py = pybind11;

class __attribute__((visibility("default"))) Field {
    py::object field_obj;
public:
    Field();
    Field(vector<Vec> points, vector<float> values, int dimension, bool is_complex, bool is_vector);
    Field(vector<Vec> points, vector<complex<float>> values, int dimension = 3, bool is_complex = false, bool is_vector = false);
    Field(const std::string &filename, int dimension = 3, bool is_complex = false, bool is_vector = false);
    float operator()(std::vector<double> coords);
    float operator()(Vec coords);
};

