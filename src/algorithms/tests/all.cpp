#include "src/config/load/c_config.h"
#include "linalg_tests.hpp"
#include "fft_tests.hpp"
#include "sym_tests.hpp"

extern "C" bool algorithm_tests() {
  printf("\nRunning Algorithm tests\n");

  int num_tests = 3;
  bool all_tests[num_tests] = {
      linalg_tests(),
      fft_tests(),
      sym_tests(),
  };
  return print_test_results(all_tests, num_tests, "Algorithm tests");
}
