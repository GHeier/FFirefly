#include "../../config/load/c_config.h"
#include "linalg_tests.hpp"

extern "C" bool algorithm_tests() {

  int num_tests = 1;
  bool all_tests[num_tests] = {
      linalg_tests(), 
  };
  return print_test_results(all_tests, num_tests, "Algorithm tests");
}
