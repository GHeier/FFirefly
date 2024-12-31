#include <linux/limits.h>
#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <Python.h>
#include <unistd.h>

#include "field_wrapper.hpp"

namespace py = pybind11;


FieldWrapper::FieldWrapper(const std::string &filename, int dimension, bool is_complex, bool is_vector) {
    printf("filename: %s\n", filename.c_str());

    py::module sys = py::module::import("sys");
    sys.attr("path").attr("append")("/home/g/Research/bcs_diagonalization/fcode/objects");

    py::module field_module = py::module::import("field");
    printf("field_module: %s\n", field_module.attr("__name__").cast<std::string>().c_str());
    py::object load_field = field_module.attr("load_field_from_file");
    printf("load_field: %s\n", load_field.attr("__name__").cast<std::string>().c_str());

    field_obj = load_field(filename, dimension, is_complex, is_vector);
    printf("success\n");
}

double FieldWrapper::operator()(std::vector<double> coords) {
    py::object result = field_obj.attr("__call__")(*py::cast(coords));
    return result.cast<double>();
}
