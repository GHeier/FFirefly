#include "src/hamiltonian/DOS/config/load/cpp_config.hpp"
#include "src/hamiltonian/DOS/hamiltonian/band_structure.hpp"
#include "src/hamiltonian/DOS/objects/CMField/base_data.hpp"
#include "src/hamiltonian/DOS/objects/CMField/bands.hpp"
#include <fstream>

void get_band_min_max(float &emin, float &emax) {
    Bands band;
    int nx = k_mesh[0]; int ny = k_mesh[1]; int nz = k_mesh[2];
    if (dimension == 2) nz = 1;
    emin = 1000;
    emax = -1000;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                Vec kvec(1.0 * i / nx - 0.5, 1.0 * j / ny - 0.5,
                         1.0 * k / nz - 0.5);
                kvec = brillouin_zone * kvec;
                for (int n = 1; n <= nbnd; n++) {
                    float e = band(n, kvec);
                    if (emin > e)
                        emin = e;
                    if (emax < e)
                        emax = e;
                }
            }
        }
    }
}

void run() {
    printf("Custom DOS Calculation\n");
    float emin = 0;
    float emax = 0;
    get_band_min_max(emin, emax);
    float dx = (emax - emin) / w_pts;
    vector<float> w_points;
    printf("Calculating DOS from %.5f to %.5f with %d points\n", emin, emax,
           w_pts);
    printf("Spacing is %.5f\n", dx);
    string filename = outdir + prefix + "_DOS." + filetype;
    vector<cfloat> dos_vals;
    for (float x = emin; x <= emax; x += dx) {
        int index = (x - emin) / dx;
        cout << "\rDOS Calculations: " << index + 1 << "/" << w_pts;
        vector<Vec> FS = get_FS(x);
        float DOS = get_DOS(FS);
        dos_vals.push_back(cfloat(DOS, 0));
        w_points.push_back(x);
    }

    BaseData::DataVariant data = dos_vals;
    save_data(filename, data, false, {}, {{}}, w_points);
    printf("\nSaved DOS to %s\n", filename.c_str());
}
