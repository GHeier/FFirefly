#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace py = pybind11;

class __attribute__((visibility("default"))) FieldWrapper {
    py::object field_obj;
public:
    FieldWrapper(const std::string &filename, int dimension = 3, bool is_complex = false, bool is_vector = false);
    double operator()(std::vector<double> coords);
};

