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
    call_flex2();
  }
}

void call_flex() {
  CMField vertex = load_CMField(outdir + prefix + "_chi.dat");
  float U = onsite_U;

  for (int i = 0; i < vertex.data.values.size(); i++) {
    complex<Vec> temp = vertex.data.values[i];
    complex<float> V = 0.0;
    float val_real = temp.real().x;
    float val_imag = temp.imag().x;
    if (vertex.data.is_complex) {
      complex<float> X = complex<float>(static_cast<float>(val_real),
                                        static_cast<float>(val_imag));
      if (abs(U * X) >= 1) {
        printf("Geometric series not convergent: U*X = %f\n", U * X.real());
        exit(1);
      }
      V = (U * U * X) / complex<float>(1.0f - U * X) +
          (U * U * U * X * X) / complex<float>(1.0f - U * U * X * X);
    } else {
      float X = val_real;
      if (abs(U * X) >= 1) {
        printf("Geometric series not convergent: U*X = %f\n", U * X);
        exit(1);
      }
      V = U * U * X / (1 - U * X) + U * U * U * X * X / (1 - U * U * X * X);
    }
    vertex.data.values[i] = complex<Vec>(Vec(V.real()), Vec(0.0));
  }
  save_CMField(outdir + prefix + "_vertex.dat", vertex);
  cout << "Saved to " << outdir + prefix + "_vertex.dat" << endl;
}

void call_flex2() {
    string filename = outdir + prefix + "_chi.dat";
    printf("Reading chi from %s\n", filename.c_str());
    Field_C chi(outdir + prefix + "_chi.dat");
    float U = onsite_U;
    float nx = q_mesh[0], ny = q_mesh[1], nz = q_mesh[2];
    int chidim = chi.cmf.data.dimension;
    if (chidim == 2) nz = 1;

    vector<Vec> points;
    vector<float> wpts = chi.cmf.data.w_points;
    vector<complex<Vec>> values;

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
                    values.push_back(complex<Vec>(Vec(val.real()), Vec(val.imag())));
                    if (abs(U * X) >= 1) {
                        printf("Geometric series not convergent: U*X = %f\n", U * X.real());
                        exit(1);
                    }
                }
            }
        }
    }
    printf("Saving Vertex\n");
    save_to_file(outdir + prefix + "_vertex.dat", points, values, chi.cmf.data.dimension, chi.cmf.data.with_w, chi.cmf.data.with_n, chi.cmf.data.is_complex, chi.cmf.data.is_vector);
    cout << "Saved to " << outdir + prefix + "_vertex.dat" << endl;
}
