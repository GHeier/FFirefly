#include "bands_tests.hpp"
#include "../../../config/load/c_config.h"
#include "../../vec.hpp"
#include "../bands.hpp"
#include "../cmfield.hpp"

int bpnts = 5;

float bindf(int i) { return 2.0 * i / bpnts - 1.0; }

complex<Vec> func_b(Vec x, int dim) {
  Vec val = Vec(x.x + x.y + x.z);
  val.dimension = dim;
  complex<Vec> r_val = complex<Vec>(val, val);
  return r_val;
}

bool breadwrite_2d() {
  vector<Vec> points;
  vector<complex<Vec>> values;
  for (int i = 0; i < bpnts; i++) {
    for (int j = 0; j < bpnts; j++) {
      Vec x(bindf(i), bindf(j));
      x.n = 1;
      points.push_back(x);
      values.push_back(func_b(x, 2));
    }
  }
  float test_pnt = points[2](0);
  float test_val = real(values[2])(0);

  CMData data(points, values, 2, false, true, false, false);

  float data_pnt = data.points[2](0);
  float data_val = real(data.values[2])(0);

  if (abs(test_pnt - data_pnt) > 1e-6)
    return false;
  if (abs(test_val - data_val) > 1e-6)
    return false;

  save(data, "./sample_bands.dat");
  Bands loaddata;

  Vec load_pnt = points[2];
  float load_val = loaddata(1, load_pnt);

  if (abs(test_val - load_val) > 1e-6)
    return false;

  return true;
}

bool bands_tests() {
  int num_tests = 1;
  bool all_tests[num_tests] = {
      breadwrite_2d(),
  };
  // remove("./sample_bands.dat");
  return print_test_results(all_tests, num_tests, "Bands tests");
}
