#include "../../objects/matrix.hpp"
#include "../../objects/eigenvec.hpp"
#include "../../config/load/c_config.h"
#include "../linear_algebra.hpp"
#include <vector>

#include "linalg_tests.hpp"
bool test_3_by_3() {
    Matrix A(3);
    A(0,0) = 1.0;
    A(1,1) = 2.1;
    A(2,2) = 3.2;
    Eigenvector result = power_iteration(A);

    if (abs(result.eigenvalue - 3.2) > 1e-2)
        return false;
    if (abs(result.norm() - 1.0) > 1e-4)
        return false;

    vector<Eigenvector> vecs = power_iteration(A, 0.001);

    if (abs(vecs[0].eigenvalue - 3.2) > 1e-2) 
        return false;

    Eigenvector expected(3);
    expected.eigenvector[0] = 0.0;
    expected.eigenvector[1] = 0.0;
    expected.eigenvector[2] = 1.0;

    if ((vecs[0] - expected).norm() > 1e-2)
        return false;
    return true;
}

bool linalg_tests() {
  int num_tests = 1;
  bool all_tests[num_tests] = {
      test_3_by_3(), 
  };
  return print_test_results(all_tests, num_tests, "Linear Algebra tests");
}
