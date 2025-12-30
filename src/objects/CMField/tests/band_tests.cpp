#include "band_tests.hpp"
#include "src/objects/CMField/bands.hpp"
#include "src/objects/vec.hpp"
#include "src/hamiltonian/models/band_structure.hpp"
#include <filesystem>
#include <iostream>

#include "src/config/load/c_config.h"
#include "src/config/load/cpp_config.hpp"

using namespace std;

float test_2D_ek(Vec k) {
    return -2 * t0 * (cos(k(0)) + cos(k(1)));
}
float test_3D_ek(Vec k) {
    return -2 * t0 * (cos(k(0)) + cos(k(1)) + cos(k(2)));
}

bool test_3D_TB() {
    set_global(band, "tight_binding");
    set_global(celltype, "SC");
    Bands band;
    // Dimension = 3
    Vec k0(0.0, 0.0, 0.0);
    Vec k1(0.5, 0.0, 0.0);
    Vec k2(0.0, 0.5, 0.0);
    Vec k3(0.5, 0.5, 0.0);
    Vec k4(-M_PI, -M_PI, -M_PI);

    float e0 = band(1, k0);
    float e1 = band(1, k1);
    float e2 = band(1, k2);
    float e3 = band(1, k3);
    float e4 = band(1, k4);

    return (fabs(e0 - test_3D_ek(k0)) < 1e-6 &&
            fabs(e1 - test_3D_ek(k1)) < 1e-6 &&
            fabs(e2 - test_3D_ek(k2)) < 1e-6 &&
            fabs(e3 - test_3D_ek(k3)) < 1e-6 &&
            fabs(e4 - test_3D_ek(k4)) < 1e-6);
}

bool band_tests() {
    int num_tests = 1;
    bool all_tests[num_tests] = {
        test_3D_TB(),
    };
    return print_test_results(all_tests, num_tests, "Band tests");
}
