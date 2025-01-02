#include <linux/limits.h>
#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <Python.h>
#include <unistd.h>

#include <vector>

#include "field.hpp"

namespace py = pybind11;

Field::Field() {
}

Field::Field(vector<Vec> points, vector<complex<float>> values, int dimension, bool is_complex, bool is_vector) {
    py::module sys = py::module::import("sys");
    sys.attr("path").attr("append")("/home/g/Research/fcode/fcode/objects");

    py::module field_module = py::module::import("field");
    py::object load = field_module.attr("load_field_from_data");

    vector<vector<double>> points_vec;
    for (Vec point : points) {
        vector<double> temp(dimension);
        for (int i = 0; i < dimension; i++) {
            temp[i] = point(i);
        }
        points_vec.push_back(temp);
    }
    vector<complex<double>> values_vec = vector<complex<double>>(values.begin(), values.end());
    field_obj = load(points_vec, values, dimension, is_complex, is_vector);
}

Field::Field(vector<Vec> points, vector<float> values, int dimension, bool is_complex, bool is_vector) {
    py::module sys = py::module::import("sys");
    sys.attr("path").attr("append")("/home/g/Research/fcode/fcode/objects");

    py::module field_module = py::module::import("field");
    py::object load = field_module.attr("load_field_from_data");

    vector<vector<double>> points_vec;
    for (Vec point : points) {
        vector<double> temp(dimension);
        for (int i = 0; i < dimension; i++) {
            temp[i] = point(i);
        }
        points_vec.push_back(temp);
    }
    vector<double> values_vec = vector<double>(values.begin(), values.end());
    field_obj = load(points_vec, values, dimension, is_complex, is_vector);
}

Field::Field(const std::string &filename, int dimension, bool is_complex, bool is_vector) {
    printf("filename: %s\n", filename.c_str());

    py::module sys = py::module::import("sys");
    sys.attr("path").attr("append")("/home/g/Research/fcode/fcode/objects");

    py::module field_module = py::module::import("field");
    py::object load_field = field_module.attr("load_field_from_file");

    field_obj = load_field(filename, dimension, is_complex, is_vector);
}


float Field::operator()(std::vector<double> coords) {
    py::object result = field_obj.attr("__call__")(*py::cast(coords));
    return result.cast<float>();
}

float Field::operator()(Vec coords) {
    py::object result = field_obj.attr("__call__")(*py::cast(coords));
    return result.cast<float>();
}
