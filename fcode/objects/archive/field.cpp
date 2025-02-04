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
    this->dimension = 3;
    this->is_complex = false;
    this->is_vector = false;
}

Field::Field(int dimension, bool is_complex, bool is_vector) {
    this->dimension = dimension;
    this->is_complex = is_complex;
    this->is_vector = is_vector;
}

py::module Field::init() {
    py::module sys = py::module::import("sys");
    sys.attr("path").attr("append")("/home/g/Research/fcode/fcode/objects");
    py::module field_module = py::module::import("field");
    return field_module;
}

ScalarField::ScalarField() : Field() {};
ComplexField::ComplexField() : Field() {};
VectorField::VectorField() : Field() {};
ComplexVectorField::ComplexVectorField() : Field() {};

vector<vector<double>> Vec_to_doubles(vector<Vec> points) {
    vector<vector<double>> points_vec;
    for (Vec point : points) {
        vector<double> temp(3);
        for (int i = 0; i < 3; i++) {
            temp[i] = point(i);
        }
        points_vec.push_back(temp);
    }
    return points_vec;
}

vector<double> Field::Vec_to_vector(Vec &point) {
    vector<double> temp(this->dimension);
    for (int i = 0; i < this->dimension; i++) {
        temp[i] = point(i);
    }
    return temp;
}

ScalarField::ScalarField(vector<Vec> points, vector<float> values, int dimension) : Field() {
    py::module field_module = init();
    py::object load = field_module.attr("load_field_from_data");

    vector<vector<double>> points_vec = Vec_to_doubles(points);
    this->field_obj = load(points_vec, values, dimension, false, false);
    this->dimension = dimension;
}

ScalarField::ScalarField(vector<vector<double>> points, vector<float> values, int dimension) : Field() {
    py::module field_module = init();
    py::object load = field_module.attr("load_field_from_data");
    this->field_obj = load(points, values, dimension, false, false);
    this->dimension = dimension;
}

ScalarField::ScalarField(string filename, int dimension) : Field() {
    py::module field_module = init();
    py::object load = field_module.attr("load_field_from_file");
    this->field_obj = load(filename, dimension);
    this->dimension = dimension;
}

float ScalarField::operator()(vector<double> coords) {
    py::object result = this->field_obj.attr("call_real_scalar")(*py::cast(coords));
    return result.cast<float>();
}

float ScalarField::operator()(Vec &coords) {
    vector<double> vcoords = Vec_to_vector(coords);
    py::object result = this->field_obj.attr("call_real_scalar")(*py::cast(vcoords));
    return result.cast<float>();
}

ComplexField::ComplexField(vector<Vec> points, vector<complex<float>> values, int dimension) {
    py::module field_module = init();
    py::object load = field_module.attr("load_field_from_data");

    vector<vector<double>> points_vec = Vec_to_doubles(points);
    this->field_obj = load(points_vec, values, dimension, true, false);
    this->dimension = dimension;
}

ComplexField::ComplexField(vector<vector<double>> points, vector<complex<float>> values, int dimension) {
    py::module field_module = init();
    py::object load = field_module.attr("load_field_from_data");
    this->field_obj = load(points, values, dimension, true, false);
    this->dimension = dimension;
}

ComplexField::ComplexField(string filename, int dimension) {
    py::module field_module = init();
    py::object load = field_module.attr("load_field_from_file");
    this->field_obj = load(filename, dimension, true, false);
    this->dimension = dimension;
}

complex<float> ComplexField::operator()(vector<double> coords) {
    py::object result = this->field_obj.attr("call_complex_scalar")(*py::cast(coords));
    return result.cast<complex<float>>();
}

complex<float> ComplexField::operator()(Vec &coords) {
    vector<double> vcoords = Vec_to_vector(coords);
    py::object result = this->field_obj.attr("call_complex_scalar")(*py::cast(vcoords));
    return result.cast<complex<float>>();
}

VectorField::VectorField(vector<Vec> points, vector<float> values, int dimension) {
    py::module field_module = init();
    py::object load = field_module.attr("load_field_from_data");

    vector<vector<double>> points_vec = Vec_to_doubles(points);
    this->field_obj = load(points_vec, values, dimension, false, true);
    this->dimension = dimension;
}

VectorField::VectorField(vector<vector<double>> points, vector<float> values, int dimension) {
    py::module field_module = init();
    py::object load = field_module.attr("load_field_from_data");
    this->field_obj = load(points, values, dimension, false, true);
    this->dimension = dimension;
}

VectorField::VectorField(string filename, int dimension) {
    py::module field_module = init();
    py::object load = field_module.attr("load_field_from_file");
    this->field_obj = load(filename, dimension, false, true);
    this->dimension = dimension;
}

vector<float> VectorField::operator()(vector<double> coords) {
    py::object result = this->field_obj.attr("call_real_vector")(*py::cast(coords));
    return result.cast<vector<float>>();
}

vector<float> VectorField::operator()(Vec &coords) {
    vector<double> vcoords = Vec_to_vector(coords);
    py::object result = this->field_obj.attr("call_real_vector")(*py::cast(vcoords));
    return result.cast<vector<float>>();
}

ComplexVectorField::ComplexVectorField(vector<Vec> points, vector<complex<float>> values, int dimension) {
    py::module field_module = init();
    py::object load = field_module.attr("load_field_from_data");

    vector<vector<double>> points_vec = Vec_to_doubles(points);
    this->field_obj = load(points_vec, values, dimension, true, true);
    this->dimension = dimension;
}

ComplexVectorField::ComplexVectorField(vector<vector<double>> points, vector<complex<float>> values, int dimension) {
    py::module field_module = init();
    py::object load = field_module.attr("load_field_from_data");
    this->field_obj = load(points, values, dimension, true, true);
    this->dimension = dimension;
}

ComplexVectorField::ComplexVectorField(string filename, int dimension) {
    py::module field_module = init();
    py::object load = field_module.attr("load_field_from_file");
    this->field_obj = load(filename, dimension, true, true);
    this->dimension = dimension;
}

vector<complex<float>> ComplexVectorField::operator()(vector<double> coords) {
    py::object result = this->field_obj.attr("call_complex_vector")(*py::cast(coords));
    return result.cast<vector<complex<float>>>();
}

vector<complex<float>> ComplexVectorField::operator()(Vec &coords) {
    vector<double> vcoords = Vec_to_vector(coords);
    py::object result = this->field_obj.attr("call_complex_vector")(*py::cast(vcoords));
    return result.cast<vector<complex<float>>>();
}
