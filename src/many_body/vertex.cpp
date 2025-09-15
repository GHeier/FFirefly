#include <complex>
#include <vector>

#include "../config/load/cpp_config.hpp"
#include "../objects/CMField/cmfield.hpp"
#include "../objects/CMField/fields.hpp"
#include "../objects/vec.hpp"
#include "vertex.hpp"

using namespace std;

extern "C" void vertex_wrapper() {
  if (interaction == "FLEX") {
    call_flex();
  }
  else 
      printf("Interaction '%s' not available\n", interaction.c_str());
}

void call_flex() {
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
    vector<vector<vector<float>>> vec_values(1);

    printf("Computing vertex\n");
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
                    complex<float> val = (U * U * X) / complex<float>(1.0f - U * X) + (U * U * U * X * X) / complex<float>(1.0f - U * U * X * X);
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
    }
    printf("Saving Vertex\n");
    string file = outdir + prefix + "_vertex." + filetype;
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
        save_to_field(file, vec_values, BZ, mesh, chi.cmf.data.w_points, chi.cmf.data.is_complex, chi.cmf.data.is_vector);
    }
    cout << "Saved to " << outdir + prefix + "_vertex." + filetype << endl;
}
