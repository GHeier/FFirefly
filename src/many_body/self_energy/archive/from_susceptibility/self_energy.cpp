#include <complex>
#include <vector>

#include "src/config/load/cpp_config.hpp"
#include "src/objects/CMField/fields.hpp"
#include "src/objects/vec.hpp"
#include "self_energy.hpp"

using namespace std;

void call_flex_self_energy() {
    string filename = outdir + prefix + "_chi." + filetype;
    printf("Reading chi from %s\n", filename.c_str());
    Field_R chi(filename);
    float U = U0;
    int chidim = chi.cmf.data.dimension;

    vector<float> wpts = chi.cmf.data.w_points;
    if (wpts.size() == 0) {
        wpts.push_back(0.0);
    }
    vector<float> vals;

    printf("Computing vertex\n");

    float ave_slope = 0.0;

    if (!chi.cmf.data.as_mesh) {
        printf("Using stored point data\n");
        // Loop over stored points
        for (const auto& point : chi.cmf.data.points) {
            Vec q(point.data(), chidim);
            q.dimension = chidim;
            for (int l = 0; l < wpts.size(); l++) {
                float w = wpts[l];

                float X = chi(q, w);
                //cfloat val = (U * U * X) / cfloat(1.0f - U * X) + (U * U * U * X * X) / cfloat(1.0f - U * U * X * X);
                float val = (U * U * U * X * X) / (1.0f - U * X) + (U * U * X) / (1.0f - U * U * X * X);
                vals.push_back(val);

                if (l > 0 && wpts[l - 1] * wpts[l] < 0)
                    ave_slope -= (vals[vals.size() - 1] - vals[vals.size() - 2]) / (wpts[l] - wpts[l-1]);

                if (abs(U * X) >= 1) {
                    printf("Geometric series not convergent: U*X = %f\n", U * X);
                    exit(1);
                }
            }
        }
        ave_slope /= chi.cmf.data.points.size();
    } else {
        printf("Using mesh data\n");
        // Loop over mesh
        float nx = q_mesh[0], ny = q_mesh[1], nz = q_mesh[2];
        if (chidim == 2) nz = 1;

        int nk = static_cast<int>(nx * ny * nz);

        // w loop outermost to generate data in w-k order (matching CMF_search expectations)
        for (int l = 0; l < wpts.size(); l++) {
            float w = wpts[l];
            float tmp_slope = 0.0;
            for (int i = 0; i < nx; i++) {
                for (int j = 0; j < ny; j++) {
                    for (int k = 0; k < nz; k++) {
                        Vec q = brillouin_zone * Vec(i / nx - 0.5, j / ny - 0.5, k / nz - 0.5);
                        q.dimension = chidim;

                        float X = chi(q, w);
                        //cfloat val = (U * U * X) / cfloat(1.0f - U * X) + (U * U * U * X * X) / cfloat(1.0f - U * U * X * X);
                        float val = (U * U * U * X * X) / (1.0f - U * X) + (U * U * X) / (1.0f - U * U * X * X);
                        vals.push_back(val);

                        if (l > 0 && wpts[l - 1] * wpts[l] < 0) {
                            int k_idx = i * static_cast<int>(ny * nz) + j * static_cast<int>(nz) + k;
                            int prev_w_idx = (l - 1) * nk + k_idx;
                            ave_slope -= (vals[vals.size() - 1] - vals[prev_w_idx]) / (wpts[l] - wpts[l-1]);
                        }

                        if ((U * X) >= 1.0) {
                            printf("Geometric series not convergent: U*X = %f\n", U * X);
                            exit(1);
                        }
                    }
                }
            }
        }
        ave_slope /= (nx * ny * nz);
    }

    float max_val = 0.0f;
    for (const auto& v : vals) {
        if (abs(v) > max_val) {
            max_val = abs(v);
        }
    }
    printf("Computed %zu vertex values\n", vals.size());
    printf("Max Self-Energy magnitude: %f\n", max_val);
    printf("Renormalization: %f\n", 1 + ave_slope);
    printf("Quasiparticle Weight: %f\n", 1/(1 + ave_slope));

    printf("Saving Self-Energy\n");
    string file = outdir + prefix + "_self_energy." + filetype;
    if (filetype == "hdf5" || filetype == "h5") {
        if (chi.cmf.data.as_mesh) {
            vector<int> mesh_vec(q_mesh.begin(), q_mesh.begin() + chidim);
            // save_data(filename, data, inds, mesh, domain, w_points, points)
            save_data(file, vals, vector<int>{}, mesh_vec, brillouin_zone, wpts);
        } else {
            save_data(file, vals, vector<int>{}, vector<int>{}, vector<vector<float>>{{}}, wpts, chi.cmf.data.points);
        }
    }
    cout << "Saved to " << outdir + prefix + "_self_energy." + filetype << endl;
    Field_R chi2(outdir + prefix + "_chi.h5");
    Field_R sigma2(outdir + prefix + "_self_energy.h5");
    Vec q(0.6, 0.8, -0.2);
    float val_chi = (chi2(q));
    float val_sigma = (sigma2(q));

    //cfloat val = (U * U * X) / cfloat(1.0f - U * X) + (U * U * U * X * X) / cfloat(1.0f - U * U * X * X);
    float val = (U * U * U * val_chi * val_chi) / (1.0f - U * val_chi) + (U * U * val_chi) / (1.0f - U * U * val_chi * val_chi);

    printf("At Vec q = (0.6, 0.8, -0.2):\n");
    printf("Chi(q) = %f\n", val_chi);
    printf("Sigma(q) = %f\n", val_sigma);
    printf("Expected Sigma(q) = %f\n", val);
}

