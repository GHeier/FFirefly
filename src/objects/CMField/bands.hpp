#include "fields.hpp"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

class Bands {
public:
  Field_R fields;
  bool file_found;

  Bands();
  //~Bands();
  float operator()(int n, Vec k);
  float operator()(Vec k);
};

Vec vk(int n, Vec k, Bands &band);
