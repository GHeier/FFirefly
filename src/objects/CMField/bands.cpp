#include "bands.hpp"
#include "fields.hpp"
#include "../../hamiltonian/band_structure.hpp"
#include "../../config/load/cpp_config.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>


using namespace std;

Bands::Bands() {
    file_found = false;
    ifstream file(outdir + prefix + "_bands.dat");
    if (file.is_open()) {
        Field_R data(outdir + prefix + "_bands.dat");
    }
    file_found = true;
}

float Bands::operator()(int n, Vec k) {
    if (n < 1 || n > fields.size()) {
        printf("Error: band index %d out of range\n", n);
        exit(1);
    }
    n--;
    if (file_found) 
        return fields[n](k);
    return epsilon(n, k);
}

float Bands::operator()(Vec k) {
    return operator()(1, k);
}
