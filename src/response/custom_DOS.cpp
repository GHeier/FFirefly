#include "../config/load/cpp_config.hpp"
#include "../hamiltonian/band_structure.hpp"
#include "../objects/CMField/cmfield.hpp"
#include <fstream>

void get_band_min_max(float &emin, float &emax) {
    int z_num = 200;
    if (dimension == 2) {
        z_num = 1;
    }
    emin = 1000;
    emax = -1000;
    for (int i = 0; i < 200; i++) {
        for (int j = 0; j < 200; j++) {
            for (int k = 0; k < z_num; k++) {
                Vec kvec(1.0 * i / 200 - 0.5, 1.0 * j / 200 - 0.5,
                         1.0 * k / 200 - 0.5);
                kvec = brillouin_zone * kvec;
                for (int n = 1; n <= nbnd; n++) {
                    float e = epsilon(n, kvec);
                    if (emin > e)
                        emin = e;
                    if (emax < e)
                        emax = e;
                }
            }
        }
    }
}

extern "C" void surface_sum() {
    printf("Custom DOS Calculation\n");
    float emin = 0;
    float emax = 0;
    get_band_min_max(emin, emax);
    float dx = (emax - emin) / w_pts;
    vector<float> dos_list;
    vector<float> w_points;
    printf("Calculating DOS from %.5f to %.5f with %d points\n", emin, emax,
           w_pts);
    printf("Spacing is %.5f\n", dx);
    string filename = outdir + prefix + "_DOS." + filetype;
    if (filetype == "dat" || filetype == "txt") {
        ofstream file(filename);
        file << "# w f" << endl;
        for (float x = emin; x <= emax; x += dx) {
            vector<Vec> FS = get_FS(x);
            float DOS = get_DOS(FS);
            dos_list.push_back(DOS);
            file << x << " " << DOS << endl;
            int index = (x - emin) / dx;
            cout << "\rDOS Calculations: " << index + 1 << "/" << w_pts;
        }
    }
    else if (filetype == "hdf5" || filetype == "h5") {
        vector<vector<vector<float>>> dos_vals(1);
        for (float x = emin; x <= emax; x += dx) {
            vector<Vec> FS = get_FS(x);
            float DOS = get_DOS(FS);
            int index = (x - emin) / dx;
            cout << "\rDOS Calculations: " << index + 1 << "/" << w_pts;
            vector<float> dos_vec = {DOS};
            dos_vals[0].push_back(dos_vec);
            w_points.push_back(x);
        }
        vector<vector<float>> earr = {{emax - emin}};
        save_to_field(filename, dos_vals, earr, {w_pts}, w_points, false, false);
    }
    printf("\nSaved DOS to %s\n", filename.c_str());
}
