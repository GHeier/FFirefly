#include "src/config/load/cpp_config.hpp"
#include "src/hamiltonian/models/band_structure.hpp"
#include "main.hpp"
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

float self_energy_renormalization() {
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
    vector<int> mesh_for_save = sigma.cmf.data.as_mesh ?
        vector<int>{sigma.cmf.data.mesh[0], sigma.cmf.data.mesh[1], sigma.cmf.data.mesh[2]} :
        q_mesh;
    vector<vector<float>> domain_to_use = sigma.cmf.data.domain.empty() ? brillouin_zone : sigma.cmf.data.domain;
    save_data(file, vals, sigma.cmf.data.inds, mesh_for_save, domain_to_use);
    cout << "Saved to " << file << endl;

    vector<Vec> FS = get_FS(fermi_energy);
    ave = 0;
    double norm = 0;
    for (Vec k : FS) {
        double dk = vp(k.n, k) * k.area;
        ave += -imag(sigma(k, 1e-4) - sigma(k, -1e-4)) / (2e-4) * dk;
        norm += dk;
    }
    float m_star = 1 + ave / norm;
    printf("Average m* on Fermi Surface is: %lf\n", m_star);
    return m_star;
}


