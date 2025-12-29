#pragma once

#include <complex>

#include "../vec.hpp"
#include "../eigenvec.hpp"
#include "field.hpp"

using namespace std;
using cfloat = complex<float>;

// Struct for eigenvalue/eigenvector pairs
struct eigvec {
    float eigenvalue;
    vector<cfloat> eigenvector;

    eigvec(int size) : eigenvector(size) {}
};

// Scalar fields
class Field_C {
public:
  FieldImpl cmf;

  Field_C();
  Field_C(const BaseData::DataVariant& data,
          const vector<int>& mesh = {},
          const vector<vector<float>>& domain = {},
          const vector<float>& w_points = {},
          bool centered = true);
  Field_C(FieldImpl f);
  Field_C(const string& filename, bool centered = true);

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
  FieldImpl cmf;

  Field_R();
  Field_R(const BaseData::DataVariant& data,
          const vector<int>& mesh = {},
          const vector<vector<float>>& domain = {},
          const vector<float>& w_points = {},
          bool centered = true);
  Field_R(FieldImpl f);
  Field_R(const string& filename, bool centered = true);

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
  FieldImpl cmf;

  Field_CM();
  Field_CM(const BaseData::DataVariant& data,
           const vector<int>& inds,
           const vector<int>& mesh = {},
           const vector<vector<float>>& domain = {},
           const vector<float>& w_points = {},
           bool centered = true);
  Field_CM(FieldImpl f);
  Field_CM(const string& filename, bool centered = true);

  // Copy assignment operator
  Field_CM& operator=(const Field_CM& other);

  // Save to file
  void save(const string& filename);

  // Returns full matrix at point
  vector<vector<cfloat>> operator()(Vec point, float w = 0);

  // Returns full matrix at frequency only (k-independent, averaged over k if needed)
  vector<vector<cfloat>> operator()(float w);

  // List-based operator for multiple points
  vector<vector<vector<cfloat>>> operator()(const vector<Vec>& points, float w = 0);

  // List-based operator for multiple w-points
  vector<vector<vector<cfloat>>> operator()(const vector<float>& w_points);

  // Diagonalize matrix at point and return eigenvalues
  vector<float> diag(Vec point, float w = 0);

  // Diagonalize matrix at point and return eigenvalues and eigenvectors
  vector<eigvec> fulldiag(Vec point, float w = 0);

  // Get underlying data
  BaseData* get_data();
};

class Field_RM {
public:
  FieldImpl cmf;

  Field_RM();
  Field_RM(const BaseData::DataVariant& data,
           const vector<int>& inds,
           const vector<int>& mesh = {},
           const vector<vector<float>>& domain = {},
           const vector<float>& w_points = {},
           bool centered = true);
  Field_RM(FieldImpl f);
  Field_RM(const string& filename, bool centered = true);

  // Copy assignment operator
  Field_RM& operator=(const Field_RM& other);

  // Save to file
  void save(const string& filename);

  // Returns full matrix at point
  vector<vector<float>> operator()(Vec point, float w = 0);

  // Returns full matrix at frequency only (k-independent, averaged over k if needed)
  vector<vector<float>> operator()(float w);

  // List-based operator for multiple points
  vector<vector<vector<float>>> operator()(const vector<Vec>& points, float w = 0);

  // List-based operator for multiple w-points
  vector<vector<vector<float>>> operator()(const vector<float>& w_points);

  // Diagonalize matrix at point and return eigenvalues
  vector<float> diag(Vec point, float w = 0);

  // Diagonalize matrix at point and return eigenvalues and eigenvectors
  vector<eigvec> fulldiag(Vec point, float w = 0);

  // Get underlying data
  BaseData* get_data();
};

// Unified Field with runtime type dispatch
class Field {
public:
  bool is_complex;
  bool is_vector;
  bool is_matrix;

  Field_C* field_c;
  Field_R* field_r;
  Field_CM* field_cm;
  Field_RM* field_rm;

  // Plot metadata
  string default_plot_type;
  string title;
  string x_label;
  string y_label;

  Field();
  Field(const string& filename, bool centered = true);
  ~Field();

  // Generate plot labels from filename
  void generate_plot_labels(const string& filename);

  // Save to file
  void save(const string& filename);

  // Scalar complex operator
  complex<float> operator_scalar_complex(Vec point, float w = 0);
  complex<float> operator_scalar_complex(float w);
  vector<complex<float>> operator_scalar_complex(const vector<Vec>& points, float w = 0);
  vector<complex<float>> operator_scalar_complex(const vector<float>& w_points);

  // Scalar real operator
  float operator_scalar_real(Vec point, float w = 0);
  float operator_scalar_real(float w);
  vector<float> operator_scalar_real(const vector<Vec>& points, float w = 0);
  vector<float> operator_scalar_real(const vector<float>& w_points);

  // Matrix complex operator
  vector<vector<complex<float>>> operator_matrix_complex(Vec point, float w = 0);
  vector<vector<vector<complex<float>>>> operator_matrix_complex(const vector<Vec>& points, float w = 0);

  // Matrix real operator
  vector<vector<float>> operator_matrix_real(Vec point, float w = 0);
  vector<vector<vector<float>>> operator_matrix_real(const vector<Vec>& points, float w = 0);

  // Get underlying data
  BaseData* get_data();
};
vector<eigvec> fulldiag(vector<vector<cfloat>> &matrix);
vector<float> diag(vector<vector<cfloat>> &matrix);
