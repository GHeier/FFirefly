#include <string>
#include <vector>

#include "../../config/load/c_config.h"
#include "../CMData/tests/all.hpp"
//#include "../CMField/tests/bands_tests.hpp"
//#include "../CMField/tests/cmfield_tests.hpp"
#include "../CMField/tests/base_data_tests.hpp"
#include "../CMField/tests/field_tests.hpp"
#include "../CMField/tests/field_wrapper_tests.hpp"
#include "../CMField/tests/matrix_field_tests.hpp"
#include "../CMField/tests/band_tests.hpp"
#include "all.hpp"
#include "surface_tests.hpp"
#include "array_tests.hpp"

using namespace std;

extern "C" bool object_tests() {
    printf("\nRunning Object tests\n");
    int num_tests = 8;
    bool all_tests[num_tests] = {
        base_data_tests(),
        field_tests(),
        field_wrapper_tests(),
        matrix_field_tests(),
        CMData_tests(),
        //cmfield_tests(),
        band_tests(),
        surface_tests(),
        array_tests(),
    };
    return print_test_results(all_tests, num_tests, "Object tests");
}
