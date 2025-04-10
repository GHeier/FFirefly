#include "bands.hpp"
#include "../../config/load/cpp_config.hpp"
#include "../../hamiltonian/band_structure.hpp"
#include "fields.hpp"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

Bands::Bands() {
  file_found = false;
  string filename = outdir + prefix + "_bands.dat";
  ifstream file(filename);
  if (file.is_open()) {
    Field_R data(filename);
    fields = data;
    file_found = true;
    printv("Read in bands from %s\n", filename);
  } else
    printv("Taking bands from .cfg file; No %s band file found\n", filename);
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
