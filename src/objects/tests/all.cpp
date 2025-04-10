#include <string>
#include <vector>

#include "../../config/load/c_config.h"
#include "../../config/load/jl_interface.h"
#include "../CMData/tests/all.hpp"
#include "../CMField/tests/bands_tests.hpp"
#include "../CMField/tests/cmfield_tests.hpp"
#include "all.hpp"

using namespace std;

extern "C" bool object_tests() {
  printf("\nRunning Object tests\n");
  int num_tests = 3;
  bool all_tests[num_tests] = {
      // field_tests(),
      CMData_tests(),
      cmfield_tests(),
      bands_tests(),
  };
  return print_test_results(all_tests, num_tests, "Object tests");
}

bool field_tests() {
  string folder = "objects/tests/";
  string file = "field_tests";
  string module = "FieldTests";

  // string function = "test_2dim_real_scalar_field";
  // call_julia_func(folder.c_str(), file.c_str(), module.c_str(),
  // function.c_str());
  return false;
}
