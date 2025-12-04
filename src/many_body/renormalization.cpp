#include "../config/load/cpp_config.hpp"
#include "../hamiltonian/band_structure.hpp"
#include "renormalization.hpp"
#include "../objects/CMField/fields.hpp"

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
    string filename = outdir + prefix + "_sigma." + filetype;
    printf("Reading self_energy from %s\n", filename.c_str());
    Field_C sigma(filename);

    vector<cfloat> vals;
    float maxval = 0;
    double ave = 0.0;

    if (!sigma.cmf.data.as_mesh) {
        // Loop over stored points
        for (const auto& point : sigma.cmf.data.points) {
            Vec kvec(point.data(), dimension);
            kvec.dimension = dimension;
            float slope_r = real(sigma(kvec, 1e-4) - sigma(kvec, -1e-4)) / (2e-4);
            float slope_i = imag(sigma(kvec, 1e-4) - sigma(kvec, -1e-4)) / (2e-4);
            float slope = slope_r + slope_i;
            if (maxval < fabs(slope))
                maxval = fabs(slope);
            vals.push_back(cfloat(-slope_i, 0.0f));
            ave += fabs(slope_i) / sigma.cmf.data.points.size();
        }
    } else {
        // Loop over mesh
        int kx = k_mesh[0];
        int ky = k_mesh[1];
        int kz = k_mesh[2];
        if (dimension == 2)
            kz = 1;

        for (int i = 0; i < kx; i++) {
            for (int j = 0; j < ky; j++) {
                for (int k = 0; k < kz; k++) {
                    Vec kvec = get_kvec(i, j, k);
                    float slope_r = real(sigma(kvec, 1e-4) - sigma(kvec, -1e-4)) / (2e-4);
                    float slope_i = imag(sigma(kvec, 1e-4) - sigma(kvec, -1e-4)) / (2e-4);
                    float slope = slope_r + slope_i;
                    if (maxval < fabs(slope))
                        maxval = fabs(slope);
                    vals.push_back(cfloat(-slope_i, 0.0f));
                    ave += fabs(slope_i) / (kx * ky * kz);
                    if (dimension == 2)
                        break;
                }
            }
        }
    }

    printf("Max m*(q) = %f\n", 1 + maxval);
    printf("Average m* on sampled points: %f\n", 1 + ave);
    printf("Saving Renormalization\n");
    string file = outdir + prefix + "_renormalization." + filetype;
    BaseData::DataVariant dv = vals;
    vector<int> mesh_for_save = sigma.cmf.data.as_mesh ?
        vector<int>{sigma.cmf.data.mesh[0], sigma.cmf.data.mesh[1], sigma.cmf.data.mesh[2]} :
        q_mesh;
    vector<vector<float>> domain_to_use = sigma.cmf.data.domain.empty() ? brillouin_zone : sigma.cmf.data.domain;
    save_data(file, dv, true, mesh_for_save, domain_to_use);
    cout << "Saved to " << file << endl;

    vector<Vec> FS = get_FS(fermi_energy);
    ave = 0;
    double norm = 0;
    for (Vec k : FS) {
        double dk = vp(k.n, k) * k.area;
        ave += -imag(sigma(k, 1e-4) - sigma(k, -1e-4)) / (2e-4) * dk;
        norm += dk;
    }
    printf("Average m* on Fermi Surface is: %lf\n", 1 + ave / norm);
}

void FLEX_renormalization() {
    string filename = outdir + prefix + "_chi." + filetype;
    printf("Reading chi from %s\n", filename.c_str());
    Field_C chi(filename);
    float U = onsite_U;
    int chidim = chi.cmf.data.dimension;

    vector<cfloat> vals;
    float maxval = 0;
    float ave = 0.0;

    printf("Computing Renormalization\n");

    if (!chi.cmf.data.as_mesh) {
        printf("Using stored points for renormalization calculation\n");
        // Loop over stored points
        for (const auto& point : chi.cmf.data.points) {
            Vec q(point.data(), chidim);
            q.dimension = chidim;

            cfloat X = chi(q);
            cfloat val = (U*U*U * X*X) / cfloat(1.0f - U * X) + (U*U * X) / cfloat(1.0f - U * U * X * X);
            vals.push_back(val);
            ave += val.real() / chi.cmf.data.points.size();

            if (maxval < val.real())
                maxval = val.real();
            if (abs(U * X) >= 1) {
                printf("Geometric series not convergent: U*X = %f\n", U * X.real());
                exit(1);
            }
        }
    } else {
        printf("Using mesh for renormalization calculation\n");
        // Loop over mesh - use chi's mesh if available, otherwise use q_mesh
        vector<int> mesh_to_use = chi.cmf.data.as_mesh ?
            vector<int>{chi.cmf.data.mesh[0], chi.cmf.data.mesh[1], chi.cmf.data.mesh[2]} :
            q_mesh;

        float nx = mesh_to_use[0], ny = mesh_to_use[1], nz = mesh_to_use[2];
        if (chidim == 2) nz = 1;

        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                for (int k = 0; k < nz; k++) {
                    Vec q = brillouin_zone * Vec(i / nx - 0.5, j / ny - 0.5, k / nz - 0.5);
                    q.dimension = chidim;

                    cfloat X = chi(q);
                    cfloat val = (U*U*U * X*X) / cfloat(1.0f - U * X) + (U*U * X) / cfloat(1.0f - U * U * X * X);
                    vals.push_back(val);
                    //vals.push_back(1.0);
                    ave += val.real() / (nx * ny * nz);

                    if (maxval < val.real())
                        maxval = val.real();
                    if (abs(U * X) >= 1) {
                        printf("Geometric series not convergent: U*X = %f\n", U * X.real());
                        exit(1);
                    }
                }
            }
        }
    }

    printf("Max m*(q) = %f\n", 1 + maxval);
    printf("Average m* on sampled points: %f\n", 1 + ave);
    printf("Saving Renormalization\n");
    string file = outdir + prefix + "_renormalization." + filetype;

    // Use chi's domain and mesh to ensure consistency
    BaseData::DataVariant dv = vals;
    vector<int> mesh_for_save = chi.cmf.data.as_mesh ?
        vector<int>{chi.cmf.data.mesh[0], chi.cmf.data.mesh[1], chi.cmf.data.mesh[2]} :
        q_mesh;
    vector<vector<float>> domain_to_use = chi.cmf.data.domain.empty() ? brillouin_zone : chi.cmf.data.domain;
    save_data(file, dv, true, mesh_for_save, domain_to_use);
    //if (filetype == "hdf5" || filetype == "h5") {
    //    chi.cmf.data.data = vals;
    //    chi.save(file);
    //}
    cout << "Saved to " << outdir + prefix + "_renormalization." + filetype << endl;

}


