#include "../config/load/cpp_config.hpp"
#include "../hamiltonian/band_structure.hpp"
#include "renormalization.hpp"
#include "../objects/CMField/fields.hpp"

extern "C" void renormalization_wrapper() {
    if (method == "analytic") {
        if (interaction == "FLEX")
            FLEX_renormalization();
        else
            printf("Analytic renormalization for '%s' interaction not available\n", interaction.c_str());
    }
    else 
        self_energy_renormalization();
}

Vec get_kvec(int i, int j, int k) {
    Vec v(
        1.0 * i / k_mesh[0] - 0.5, 
        1.0 *j / k_mesh[1] - 0.5, 
        1.0 * k / k_mesh[2] - 0.5,
        3
    );
    Vec newv = brillouin_zone * v;
    return newv;
}

void self_energy_renormalization() {
    string filename = outdir + prefix + "_self_energy." + filetype;
    printf("Reading self_energy from %s\n", filename.c_str());
    Field_C sigma(filename);
    int kx = k_mesh[0];
    int ky = k_mesh[1];
    int kz = k_mesh[2];
    if (dimension == 2)
        kz = 1;
    vector<vector<vector<float>>> eff_mass(1);
    float maxval = 0;
    for (int i = 0; i < kx; i++) {
        for (int j = 0; j < ky; j++) {
            for (int k = 0; k < kz; k++) {
                Vec kvec = get_kvec(i, j, k);
                float slope_r = real(sigma(kvec, 1e-4) - sigma(kvec, -1e-4)) / (2e-4);
                float slope_i = imag(sigma(kvec, 1e-4) - sigma(kvec, -1e-4)) / (2e-4);
                float slope = slope_r + slope_i;
                if (maxval < fabs(slope))
                    maxval = fabs(slope);
                eff_mass[0].push_back({-slope});
                if (dimension == 2)
                    break;
            }
        }
    }
    printf("Max m*(q) = %f\n", 1 + maxval);
    printf("Saving Renormalization\n");
    string file = outdir + prefix + "_renormalization." + filetype;
    //if (filetype == "dat" || filetype == "txt")
    //    save_to_file(file, points, values, chi.cmf.data.dimension, chi.cmf.data.with_w, chi.cmf.data.with_n, chi.cmf.data.is_complex, chi.cmf.data.is_vector);
    if (filetype == "hdf5" || filetype == "h5") {
        vector<vector<float>> BZ = brillouin_zone;
        BZ.resize(dimension);
        for (auto &row : BZ)
            row.resize(dimension);
        Vec first = BZ * Vec(-0.5, -0.5, -0.5);
        vector<int> mesh = k_mesh;
        if (dimension == 2)
            mesh = {k_mesh[0], k_mesh[1]};
        save_to_field(file, eff_mass, BZ, mesh, {}, false, false);
    }
    cout << "Saved to " << outdir + prefix + "_renormalization." + filetype << endl;

    vector<Vec> FS = get_FS(fermi_energy);
    double ave = 0;
    double norm = 0;
    for (Vec k : FS) {
        double dk = vp(k.n, k) * k.area;
        ave += -real(sigma(k, 1e-4) - sigma(k, -1e-4)) / (2e-4) * dk;
        norm += dk;
    }
    printf("Average m* on Fermi Surface is: %lf\n", 1 + ave / norm);
}

void FLEX_renormalization() {
    string filename = outdir + prefix + "_chi." + filetype;
    printf("Reading chi from %s\n", filename.c_str());
    Field_C chi(filename);
    float U = onsite_U;
    float nx = q_mesh[0], ny = q_mesh[1], nz = q_mesh[2];
    int chidim = chi.cmf.data.dimension;
    if (chidim == 2) nz = 1;

    vector<Vec> points;
    vector<complex<Vec>> values;
    vector<vector<vector<float>>> vec_values(1);

    float maxval = 0;
    printf("Computing Renormalization\n");
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                Vec q = brillouin_zone * Vec(i / nx - 0.5, j / ny - 0.5, k / nz - 0.5);
                q.dimension = chidim;
                points.push_back(q);
                complex<float> X = chi(q);
                complex<float> val = (U*U*U * X*X) / complex<float>(1.0f - U * X) + (U*U * X) / complex<float>(1.0f - U * U * X * X);
                if (maxval < val.real())
                    maxval = val.real();
                if (filetype == "dat" || filetype == "txt")
                    values.push_back(complex<Vec>(Vec(val.real()), Vec(val.imag())));
                else if (filetype == "h5" || filetype == "hdf5")
                    vec_values[0].push_back({val.real(), val.imag()});
                if (abs(U * X) >= 1) {
                    printf("Geometric series not convergent: U*X = %f\n", U * X.real());
                    exit(1);
                }
            }
        }
    }
    printf("Max m*(q) = %f\n", 1 + maxval);
    printf("Saving Renormalization\n");
    string file = outdir + prefix + "_renormalization." + filetype;
    if (filetype == "dat" || filetype == "txt")
        save_to_file(file, points, values, chi.cmf.data.dimension, chi.cmf.data.with_w, chi.cmf.data.with_n, chi.cmf.data.is_complex, chi.cmf.data.is_vector);
    else if (filetype == "hdf5" || filetype == "h5") {
        vector<vector<float>> BZ = brillouin_zone;
        BZ.resize(dimension);
        for (auto &row : BZ)
            row.resize(dimension);
        Vec first = BZ * Vec(-0.5, -0.5, -0.5);
        vector<int> mesh = q_mesh;
        if (chi.cmf.data.dimension == 2)
            mesh = {q_mesh[0], q_mesh[1]};
        save_to_field(file, vec_values, BZ, mesh, {}, false, false);
    }
    cout << "Saved to " << outdir + prefix + "_renormalization." + filetype << endl;

}


