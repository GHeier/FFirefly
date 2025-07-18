#include "bands.hpp"
#include "../../config/load/cpp_config.hpp"
#include "../../hamiltonian/band_structure.hpp"
#include "../vec.hpp"
#include "fields.hpp"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

Bands::Bands() {
    file_found = false;
    string filename = outdir + prefix + "_bands." + filetype;
    ifstream file(filename);
    if (file.is_open() and automatic_file_read) {
        Field_R data(filename);
        fields = data;
        file_found = true;
        printv("Read in bands from %s\n", filename);
    } else
        printv("Taking bands from .cfg file; No %s band file found\n",
               filename.c_str());
}

float Bands::operator()(int n, Vec k) {
    if (n < 1 || (file_found and n > fields.cmf.values.size())) {
        printf("Error: band index %d out of range\n", n);
        exit(1);
    }
    if (file_found)
        return fields(n, k, 0.0);
    return epsilon(n, k);
}

float Bands::operator()(Vec k) { return operator()(1, k); }


Vec vk(int n, Vec k, Bands &band) {
    Vec dk = brillouin_zone * Vec(1e-3, 1e-3, 1e-3);
    Vec dx = Vec(dk.x);
    Vec dy = Vec(0, dk.y);
    Vec dz = Vec(0, 0, dk.z);
    float dfdx = (band(n, k + dx) - band(n, k - dx)) / (2*dk.x);
    float dfdy = (band(n, k + dy) - band(n, k - dy)) / (2*dk.y);
    float dfdz = (band(n, k + dz) - band(n, k - dz)) / (2*dk.z);
    return Vec(dfdx, dfdy, dfdz);
}
