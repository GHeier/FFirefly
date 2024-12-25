#include <math.h>
#include <vector>
#include "../config/load/cpp_config.hpp"

int n = k_mesh[0];
int m = q_mesh[0];
int l = 5;
float wc = bcs_cutoff_frequency;
float mu = fermi_energy;
float U = onsite_U;
// Constants
float t = 1.0;
float tn = 0.0;
float tnn = 0.0;
float k_max = M_PI;
int dim = dimension;

void load_cpp_cfg() {
    n = k_mesh[0];
    m = q_mesh[0];
    l = 5;
    wc = bcs_cutoff_frequency;
    mu = fermi_energy;
    U = onsite_U;
    // Constants
    t = 1.0;
    tn = 0.0;
    tnn = 0.0;
    k_max = M_PI;
    dim = dimension;
}

