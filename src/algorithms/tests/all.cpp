#include "../../config/load/c_config.h"
#include "linalg_tests.hpp"
#include "fft_tests.hpp"

extern "C" bool algorithm_tests() {
  printf("\nRunning Algorithm tests\n");

  int num_tests = 2;
  bool all_tests[num_tests] = {
      linalg_tests(),
      fft_tests(),
  };
  return print_test_results(all_tests, num_tests, "Algorithm tests");
}
