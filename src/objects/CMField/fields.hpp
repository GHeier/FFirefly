#pragma once

#include <complex>

#include "../vec.hpp"
#include "field.hpp"

using namespace std;

class Field_C {
public:
  Field cmf;

  Field_C();
  Field_C(const BaseData::DataVariant& data,
          const vector<int>& mesh = {},
          const vector<vector<float>>& domain = {},
          const vector<float>& w_points = {});
  Field_C(Field f);
  Field_C(const string& filename);

  // Copy assignment operator
  Field_C& operator=(const Field_C& other);

  // Save to file
  void save(const string& filename);

  complex<float> operator()(Vec point, float w = 0);
  complex<float> operator()(float w);
};

class Field_R {
public:
  Field cmf;

  Field_R();
  Field_R(const BaseData::DataVariant& data,
          const vector<int>& mesh = {},
          const vector<vector<float>>& domain = {},
          const vector<float>& w_points = {});
  Field_R(Field f);
  Field_R(const string& filename);

  // Copy assignment operator
  Field_R& operator=(const Field_R& other);

  // Save to file
  void save(const string& filename);

  float operator()(Vec point, float w = 0);
  float operator()(float w);
};
