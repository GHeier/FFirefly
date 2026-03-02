#include <string>
#include <vector>

#include "src/config/load/c_config.h"
#include "src/objects/CMData/tests/all.hpp"
//#include "src/objects/CMField/tests/bands_tests.hpp"
//#include "src/objects/CMField/tests/cmfield_tests.hpp"
#include "src/objects/CMField/tests/base_data_tests.hpp"
#include "src/objects/CMField/tests/field_tests.hpp"
#include "src/objects/CMField/tests/field_wrapper_tests.hpp"
#include "src/objects/CMField/tests/matrix_field_tests.hpp"
#include "src/objects/CMField/tests/hamiltonian_tests.hpp"
#include "src/objects/CMField/tests/tensor_field_tests.hpp"
#include "src/objects/CMField/tests/inds_field_tests.hpp"
#include "src/objects/CMField/tests/band_tests.hpp"
#include "src/objects/CMField/tests/quad_field_tests.hpp"
#include "all.hpp"
#include "surface_tests.hpp"
#include "array_tests.hpp"
#include "hmatrix_tests.hpp"

using namespace std;

extern "C" bool object_tests() {
    printf("\nRunning Object tests\n");
    int num_tests = 13;
    bool all_tests[num_tests] = {
        base_data_tests(),
        field_tests(),
        quad_field_tests(),
        field_wrapper_tests(),
        matrix_field_tests(),
        tensor_field_tests(),
        CMData_tests(),
        //cmfield_tests(),
        band_tests(),
        surface_tests(),
        array_tests(),
        inds_tests(),
        hamiltonian_tests(),
        hmatrix_tests(),
    };

    return print_test_results(all_tests, num_tests, "Object tests");
}
