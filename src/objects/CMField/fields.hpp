#pragma once

#include <complex>

#include "../vec.hpp"
#include "cmfield.hpp"

using namespace std;
class Field_C {
public:
  CMF cmf;

  Field_C();
  Field_C(CMF cmf);
  Field_C(string filename);

  complex<float> operator()(Vec point, float w = 0);
  complex<float> operator()(float w);
  complex<float> operator()(int n, double w);
  complex<float> operator()(int n, Vec point, float w = 0);
};

class Field_R {
public:
  CMF cmf;

  Field_R();
  Field_R(CMF cmf);
  Field_R(string filename);

  float operator()(Vec point, float w = 0);
  float operator()(double w);
  float operator()(int n, double w);
  float operator()(int n, Vec point, float w);
};
