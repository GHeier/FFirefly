#include "all.hpp"
#include "../field_wrapper.hpp"

extern "C" bool object_tests() {
    return python_field_wrapper_test();
}


bool python_field_wrapper_test() {
    FieldWrapper test_var("/home/g/Research/bcs_diagonalization/chi_mesh_static.dat");
    return true;
}
