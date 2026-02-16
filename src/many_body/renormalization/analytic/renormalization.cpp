#include "src/config/load/cpp_config.hpp"
#include "src/hamiltonian/models/band_structure.hpp"
#include "renormalization.hpp"
#include "src/objects/CMField/fields.hpp"

static Vec get_kvec(int i, int j, int k) {
    Vec v(
        1.0 * i / k_mesh[0] - 0.5, 
        1.0 *j / k_mesh[1] - 0.5, 
        1.0 * k / k_mesh[2] - 0.5,
        3
    );
    Vec newv = brillouin_zone * v;
    return newv;
}

void FLEX_renormalization() {
    string filename = outdir + prefix + "_chi." + filetype;
    printf("Reading chi from %s\n", filename.c_str());
    Field_C chi(filename);
    float U = U0;
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
    vector<int> mesh_for_save = chi.cmf.data.as_mesh ?
        vector<int>{chi.cmf.data.mesh[0], chi.cmf.data.mesh[1], chi.cmf.data.mesh[2]} :
        q_mesh;
    vector<vector<float>> domain_to_use = chi.cmf.data.domain.empty() ? brillouin_zone : chi.cmf.data.domain;
    save_data(file, vals, chi.cmf.data.inds, mesh_for_save, domain_to_use);
    //if (filetype == "hdf5" || filetype == "h5") {
    //    chi.cmf.data.data = vals;
    //    chi.save(file);
    //}
    cout << "Saved to " << outdir + prefix + "_renormalization." + filetype << endl;

}


