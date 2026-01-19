#include "run.hpp"
#include "../../../config/load/cpp_config.hpp"
#include "../../../config/load/c_config.h"
#include "../../../objects/CMField/fields.hpp"
#include <iostream>
#include <cstdio>
#include <cstdlib>

#include "../../../objects/vec.hpp"
#include "../../../hamiltonian/models/band_structure.hpp"

using cfloat = complex<float>;

void get_energy_max_min(float& E_min, float& E_max) {
    E_min = 1e10f;
    E_max = -1e10f;
    for (int i = 0; i < k_mesh[0]; i++) {
        for (int j = 0; j < k_mesh[1]; j++) {
            for (int k = 0; k < k_mesh[2]; k++) {
                Vec q = brillouin_zone * Vec(
                    (float)i / k_mesh[0] - 0.5f,
                    (float)j / k_mesh[1] - 0.5f,
                    (float)k / k_mesh[2] - 0.5f
                );
                for (int n = 1; n <= nbnd; n++) {
                    float E = epsilon(n, q);
                    if (E < E_min) E_min = E;
                    if (E > E_max) E_max = E;
                }
            }
        }
    }
    E_min += 0.01f;
    E_max -= 0.01f;
}

void electron_number(vector<cfloat> &DOS, vector<float> &w_points) {
    float dw = w_points[1] - w_points[0];
    vector<float> num_electrons(w_pts, 0.0);
    vector<cfloat> E_points(w_pts);
    for (int i = 0; i < w_pts; i++) {
        float n_e = 0.0f;
        for (int j = 0; j < i; j++) {
            n_e += 2 * real(DOS[j]) * dw;
        }
        num_electrons[i] = n_e;
        E_points[i] = cfloat(w_points[i], 0.0f);
    }
    BaseData::DataVariant dv = E_points;
    printf("Saving electron number vs mu data to %s\n", (outdir + prefix + "_n_vs_E.h5").c_str());
    save_data(outdir + prefix + "_E_vs_n.h5", dv, false, {}, {{}}, num_electrons);
}

float run() {
    vector<cfloat> DOS(w_pts, cfloat(0.0f, 0.0f));
    vector<float> w_points(w_pts);
    float E_min, E_max;
    get_energy_max_min(E_min, E_max);
    printf("E_min = %f, E_max = %f\n", E_min, E_max);
    for (int i = 0; i < w_pts; i++) {
        printf("Calculating DOS at point %d / %d\r", i + 1, w_pts);
        float E = E_min + (E_max - E_min) * i / (w_pts - 1);
        w_points[i] = E;
        vector<Vec> FS = get_FS(E);
        float sum = 0;
        for (auto& x : FS) {
            sum += x.area / vp(x.n, x);
        }
        DOS[i] = cfloat(sum / (pow(2 * M_PI, dimension)), 0.0f);
    }
    BaseData::DataVariant dv = DOS;
    printf("Saving DOS data to %s\n", (outdir + prefix + "_DOS.h5").c_str());
    save_data(outdir + prefix + "_DOS.h5", dv, false, {}, {{}}, w_points);
    electron_number(DOS, w_points);
    return real(DOS[0]);
}

int main() {
    // Load configuration from the build directory
    const char* config_path = "/home/g/Research/FFirefly/build/bin/input.cfg";

    // Load configuration using existing infrastructure
    read_c_config(config_path);
    load_cpp_config();

    // Run the method
    float result = run();

    // Success
    return 0;
}



