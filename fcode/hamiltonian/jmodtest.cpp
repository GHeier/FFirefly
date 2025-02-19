#include <string>


#include "../objects/cmf_types.hpp"  // Include CMF class
#include "../objects/vec.hpp"
#include "../config/load/c_config.h"
#include "../config/load/cpp_config.hpp"
#include "../hamiltonian/band_structure.hpp"

extern "C" double epsilon_julia(int n, double *k, int size) {
    Vec kvec;
    for (int i = 0; i < size; i++) {
        kvec(i) = k[i];
    }
    return epsilon(n, kvec);
}

extern "C" void load_config_julia(const char *filename) {
    read_c_config(filename);
    load_cpp_config();
}


extern "C" {

// Create a new CMF instance and return a pointer
CMF_CS* create_CMF_CS() {
    return new CMF_CS();
}

// Load CMF from a file
CMF_CS* load_cmf_cs(const char* filename) {
    return new CMF_CS(filename);
}

// Call operator() overload for `Vec`
void cmf_cs_call(CMF_CS* cmf, Vec point, float w, std::complex<float>* result) {
    *result = cmf->operator()(point, w);
}

// Call operator() overload for `float`
void cmf_cs_call2(CMF_CS* cmf, float w, std::complex<float>* result) {
    *result = cmf->operator()(w);
}

void cmf_save(CMF_CS* cmf, const char* filename) {
    save_CMF_to_file(filename, cmf->base);
}

// Destroy CMF instance
void destroy_CMF(CMF_CS* cmf) {
    delete cmf;
}

}
