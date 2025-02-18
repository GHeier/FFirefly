#include <iostream>
#include <string>

#include "../objects/vec.hpp"
#include "../config/load/c_config.h"
#include "../config/load/cpp_config.hpp"
#include "../hamiltonian/band_structure.hpp"

extern "C" float epsilon_julia(int n, double *k, int size) {
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
