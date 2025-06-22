#include "../config/load/cpp_config.hpp"
#include "../config/load/jl_interface.h"
#include "self_energy.hpp"
#include "../objects/CMField/fields.hpp"

extern "C" void self_energy_wrapper() {
    if (method == "sparse_ir")
        call_self_energy();
    else
        construct_self_energy();
}

void call_self_energy() {
    string folder = "hamiltonian/";
    string filename = "self_energy";
    string module = "Self_Energy";
    string function = "get_self_energy";
    call_julia_func(folder.c_str(), filename.c_str(), module.c_str(), function.c_str());
}

void construct_self_energy() {
    string filename = outdir + prefix + "_chi." + filetype;
    printf("Reading chi from %s\n", filename.c_str());
    Field_C chi(filename);
    float U = onsite_U;
    float nx = q_mesh[0], ny = q_mesh[1], nz = q_mesh[2];
    int chidim = chi.cmf.data.dimension;
    if (chidim == 2) nz = 1;

    vector<Vec> points;
    vector<float> wpts = chi.cmf.data.w_points;
    if (wpts.size() == 0) {
        wpts.push_back(0.0);
    }
    printv("wpts size: %d\n", wpts.size());
    vector<complex<Vec>> values;

    printf("Computing Self Energy\n");
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                Vec q = brillouin_zone * Vec(i / nx - 0.5, j / ny - 0.5, k / nz - 0.5);
                q.dimension = chidim;
                for (int l = 0; l < wpts.size(); l++) {
                    float w = wpts[l];
                    if (chidim == 3) q.w = w;
                    else q.z = w;
                    points.push_back(q);
                    complex<float> X = chi(q, w);
                    complex<float> val = (U*U*U * X*X) / complex<float>(1.0f - U * X) + (U*U * X) / complex<float>(1.0f - U * U * X * X);
                    values.push_back(complex<Vec>(Vec(val.real()), Vec(val.imag())));
                    if (abs(U * X) >= 1) {
                        printf("Geometric series not convergent: U*X = %f\n", U * X.real());
                        exit(1);
                    }
                }
            }
        }
    }
    printf("Saving Self Energy\n");
    save_to_file(outdir + prefix + "_self_energy." + filetype, points, values, chi.cmf.data.dimension, chi.cmf.data.with_w, chi.cmf.data.with_n, chi.cmf.data.is_complex, chi.cmf.data.is_vector);
    cout << "Saved to " << outdir + prefix + "_self_energy." + filetype << endl;
}


