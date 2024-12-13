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
