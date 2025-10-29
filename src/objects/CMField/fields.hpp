#pragma once

#include <complex>

#include "../vec.hpp"
#include "../eigenvec.hpp"
#include "field.hpp"

using namespace std;

// Scalar fields
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

  // List-based operator for multiple points
  vector<complex<float>> operator()(const vector<Vec>& points, float w = 0);

  // List-based operator for multiple w-points
  vector<complex<float>> operator()(const vector<float>& w_points);

  // Get underlying data
  BaseData* get_data();
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

  // List-based operator for multiple points
  vector<float> operator()(const vector<Vec>& points, float w = 0);

  // List-based operator for multiple w-points
  vector<float> operator()(const vector<float>& w_points);

  // Get underlying data
  BaseData* get_data();
};

// Matrix fields
class Field_CM {
public:
  Field cmf;

  Field_CM();
  Field_CM(const BaseData::DataVariant& data,
           int dim_indices,
           const vector<int>& mesh = {},
           const vector<vector<float>>& domain = {},
           const vector<float>& w_points = {});
  Field_CM(Field f);
  Field_CM(const string& filename);

  // Copy assignment operator
  Field_CM& operator=(const Field_CM& other);

  // Save to file
  void save(const string& filename);

  // Returns full matrix at point
  vector<vector<cfloat>> operator()(Vec point, float w = 0);

  // List-based operator for multiple points
  vector<vector<vector<cfloat>>> operator()(const vector<Vec>& points, float w = 0);

  // Diagonalize matrix at point and return eigenvalues
  vector<float> diag(Vec point, float w = 0);

  // Diagonalize matrix at point and return eigenvalues and eigenvectors
  vector<Eigenvector> fulldiag(Vec point, float w = 0);

  // Get underlying data
  BaseData* get_data();
};

class Field_RM {
public:
  Field cmf;

  Field_RM();
  Field_RM(const BaseData::DataVariant& data,
           int dim_indices,
           const vector<int>& mesh = {},
           const vector<vector<float>>& domain = {},
           const vector<float>& w_points = {});
  Field_RM(Field f);
  Field_RM(const string& filename);

  // Copy assignment operator
  Field_RM& operator=(const Field_RM& other);

  // Save to file
  void save(const string& filename);

  // Returns full matrix at point
  vector<vector<float>> operator()(Vec point, float w = 0);

  // List-based operator for multiple points
  vector<vector<vector<float>>> operator()(const vector<Vec>& points, float w = 0);

  // Diagonalize matrix at point and return eigenvalues
  vector<float> diag(Vec point, float w = 0);

  // Diagonalize matrix at point and return eigenvalues and eigenvectors
  vector<Eigenvector> fulldiag(Vec point, float w = 0);

  // Get underlying data
  BaseData* get_data();
};
