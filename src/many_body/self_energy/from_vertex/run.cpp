#include "run.hpp"
#include "tests/test.hpp"
#include "../../../config/load/cpp_config.hpp"
#include "../../../config/load/c_config.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "self_energy.hpp"

// Set to true to use the FFT-accelerated momentum-convolution self-energy
// calculation instead of the direct O(Nq*Nk*Nw) sum. The FFT path requires
// k_mesh == q_mesh (q_mesh is ignored, k_mesh sets both grids) and introduces
// a small Lorentzian broadening eta instead of the direct method's exact
// analytic delta-function resolution - much faster on large meshes, but an
// approximation rather than bit-identical to the direct method. See
// self_energy.cpp for the derivation of both.
bool use_fft_self_energy = false;

float run() {
    if (use_fft_self_energy)
        call_flex_self_energy_fft();
    else
        call_flex_self_energy();

    return 3.14f; // Return something of any type that can be tested in the test suite.
}

int main() {
    // Load configuration from the build directory
    std::string config_path = get_loc() + "input.cfg";

    // Load configuration using existing infrastructure
    read_c_config(config_path.c_str());
    load_cpp_config();

    // Run the method
    if (debug)
        float result = test();
    else
        float result = run();

    // Success
    return 0;
}
