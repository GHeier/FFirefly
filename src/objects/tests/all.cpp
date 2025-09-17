#include <string>
#include <vector>

#include "../../config/load/c_config.h"
#include "../CMData/tests/all.hpp"
#include "../CMField/tests/bands_tests.hpp"
#include "../CMField/tests/cmfield_tests.hpp"
#include "../CMField/tests/base_field_tests.hpp"
#include "../CMField/tests/field_tests.hpp"
#include "all.hpp"
#include "surface_tests.hpp"

using namespace std;

extern "C" bool object_tests() {
    printf("\nRunning Object tests\n");
    int num_tests = 6;
    bool all_tests[num_tests] = {
        base_field_tests(),
        CMData_tests(),
        cmfield_tests(),
        field_tests(),
        bands_tests(),
        surface_tests(),
    };
    return print_test_results(all_tests, num_tests, "Object tests");
}
