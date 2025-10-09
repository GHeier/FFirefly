#pragma once

#include <vector>
#include <string>

#include "fields.hpp"
#include "../vec.hpp"
#include "../eigenvec.hpp"

using namespace std;

class Bands {
public:
  vector<Field_R> band_fields;
  //vector<Field_RV> wavefunctions;

  bool file_found;
  int nbands;

  Bands();
  Bands(Field_CM &H);
  void fill_grid(Field_CM &H);

  // Returns eigenvalue for band n at k-point k (n starts at 1)
  float operator()(int n, Vec k);
  float operator()(Vec k);
};
