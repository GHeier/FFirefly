#include <complex>
#include <vector>

#include "../config/load/cpp_config.hpp"
#include "../objects/CondensedMatterField/fields.hpp"
#include "../objects/CondensedMatterField/cmf.hpp"
#include "../objects/vec.hpp"
#include "vertex.hpp"


using namespace std;

extern "C" void vertex_wrapper() {
    if (interaction == "FLEX") {
        call_flex();
    }
}

void call_flex() {
    CMF vertex = load_CMF_from_file(outdir + prefix + "_chi.dat"); 
    float U = onsite_U;

    for (int i = 0; i < vertex.values.size(); i++) {
        complex<Vec> temp = vertex.values[i];
        complex<float> V = 0.0;
        float val_real = temp.real().x;
        float val_imag = temp.imag().x;
        if (vertex.is_complex) {
            complex<float> X = complex<float>(static_cast<float>(val_real), static_cast<float>(val_imag));
            if (abs(U*X) >= 1) {
                printf("Geometric series not convergent: U*X = %f\n", U*X.real());
                exit(1);
            }
            V = (U*U * X) / complex<float>(1.0f - U*X) + (U*U*U * X*X) / complex<float>(1.0f - U*U * X*X);
        }
        else {
            float X = val_real;
            if (abs(U*X) >= 1) {
                printf("Geometric series not convergent: U*X = %f\n", U*X);
                exit(1);
            }
            V = U*U * X / (1 - U*X) + U*U*U * X*X / (1 - U*U * X*X);
        }
        vertex.values[i] = complex<Vec>(Vec(V.real()), Vec(0.0)); 
    }
    save_CMF_to_file(outdir+prefix+"_2PI.dat", vertex);
}
