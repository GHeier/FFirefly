#include "cmdata.hpp"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

CMData::CMData() {
  points = vector<Vec>();
  w_points = vector<float>();
  values = vector<complex<Vec>>();
  dimension = 0;
  is_complex = false;
  is_vector = false;
  with_w = false;
  with_n = false;
  n_inds = vector<int>();

  filled = false;
}

CMData::CMData(vector<Vec> points, vector<complex<Vec>> values, int dimension,
               bool with_w, bool with_n, bool is_complex, bool is_vector) {
  this->points = points;
  this->values = values;
  this->dimension = dimension;
  this->is_complex = is_complex;
  this->is_vector = is_vector;
  this->with_w = with_w;
  this->with_n = with_n;
  this->n_inds = vector<int>();

  if (with_n)
    n_inds.push_back(0);
  if (with_w or with_n) {
    for (int i = 0; i < points.size(); i++) {
      if (with_w and find(w_points.begin(), w_points.end(),
                          points[i](dimension)) == w_points.end())
        w_points.push_back(points[i](dimension));
      if (with_n and i < points.size() - 1 and
          points[i + 1].n - points[i].n != 0)
        n_inds.push_back(i + 1);
    }
  }

  filled = true;
}

CMData::CMData(string filename) { *this = load(filename); }

CMData load(string filename) {
  ifstream file(filename);
  if (!file.is_open()) {
    string base_err = "Failed to open teh file: ";
    throw runtime_error(base_err + filename);
  }
  string line;
  int ind = 0;
  int dimension = 0;
  bool w_col = false;
  bool n_col = false;
  bool is_complex = false;
  bool is_vector = false;
  vector<Vec> points;
  vector<complex<Vec>> values;
  while (getline(file, line)) {
    istringstream iss(line);
    if (ind == 0) {
      string word;
      while (iss >> word) {
        if (word == "x")
          dimension++;
        if (word == "y")
          dimension++;
        if (word == "z")
          dimension++;
        if (word == "w")
          w_col = true;
        if (word == "n")
          n_col = true;
        if (word == "Re(f)")
          is_complex = true;
        if (word == "fx")
          is_vector = true;
        if (word == "Re(fx)") {
          is_vector = true;
          is_complex = true;
        }
      }
      ind++;
      continue;
    }
    Vec point;
    float value;
    float ivalue;
    for (int i = 0; i < dimension + w_col; i++) {
      float temp;
      iss >> temp;
      point(i) = temp;
    }
    if (n_col) {
      iss >> ivalue;
      point.n = ivalue;
    }
    point.dimension = dimension;
    points.push_back(point);
    iss >> value;
    Vec f_val(value);
    Vec fi_val(value);
    complex<Vec> vval = complex<Vec>(f_val, Vec());
    if (is_complex) {
      iss >> ivalue;
      fi_val = Vec(ivalue);
      vval = complex<Vec>(f_val, fi_val);
    }
    if (is_vector) {
      f_val.dimension = dimension;
      fi_val.dimension = dimension;
      for (int i = 1; i < dimension; i++) {
        iss >> value;
        f_val(i) = value;
        if (is_complex) {
          iss >> ivalue;
          fi_val(i) = ivalue;
        }
      }
      vval = complex<Vec>(f_val, fi_val);
    }

    values.push_back(vval);
  }
  return CMData(points, values, dimension, w_col, n_col, is_complex, is_vector);
}

string get_header(int dimension, bool with_w, bool with_n, bool is_complex,
                  bool is_vector) {
  // Determine fheader
  string fheader = " f    ";
  if (is_complex) {
    fheader = "   Re(f)         Im(f)";
  }
  if (is_vector && !is_complex) {
    fheader = "   fx    ";
  }
  if (is_vector && is_complex) {
    fheader = "    Re(fx)      Im(fx) ";
  }

  // Build header
  string header = "    x        ";
  if (dimension == 0) {
    header = "    w         ";
  }
  if (dimension > 1) {
    header += "    y         ";
    if (is_vector && !is_complex) {
      fheader += "   fy        ";
    } else if (is_vector && is_complex) {
      fheader += "    Re(fy)    Im(fy) ";
    }
  }
  if (dimension > 2) {
    header += "    z        ";
    if (is_vector && !is_complex) {
      fheader += "   fz        ";
    }
    if (is_vector && is_complex) {
      fheader += "     Re(fz)    Im(fz) ";
    }
  }
  if (with_w and dimension != 0) {
    header += "   w        ";
  }
  if (with_n) {
    header += "n        ";
  }

  header += fheader;
  return header;
}

void save_to_file(string filename, vector<Vec> &points,
                  vector<complex<Vec>> &values, int dimension, bool with_w,
                  bool with_n, bool is_complex, bool is_vector) {
  ofstream file(filename);
  string header = get_header(dimension, with_w, with_n, is_complex, is_vector);
  file << fixed << setprecision(6);
  file << header << endl;
  for (int i = 0; i < points.size(); i++) {
    Vec point = points[i];
    complex<Vec> value = values[i];
    if (dimension > 0) {
      if (point(0) >= 0)
        file << " ";
      file << point(0) << "    ";
      // if (point(0) < 0) file << " ";
    }

    if (dimension > 1) {
      if (point(1) >= 0)
        file << " ";
      file << point(1) << "    ";
      // if (point(1) < 0)
      //   file << " ";
    }

    if (dimension > 2) {
      if (point(2) >= 0)
        file << " ";
      file << point(2) << "    ";
      // if (point(2) < 0)
      //   file << " ";
    }

    if (with_w) {
      if (point(dimension) >= 0)
        file << " ";
      file << point(dimension) << "    ";
      // if (point(dimension) < 0)
      //   file << " ";
    }
    if (point(dimension) >= 0)
      file << " ";
    if (with_n)
      file << point.n << "     ";

    if (is_complex) {
      file << value.real()(0) << "      " << value.imag()(0) << "      ";
      if (is_vector and dimension > 1) {
        file << value.real()(1) << "      " << value.imag()(1) << "      ";
        if (dimension > 2) {
          file << value.real()(2) << "      " << value.imag()(2) << "      ";
        }
      }
    } else {
      if (value.real()(0) >= 0)
        file << " ";
      file << value.real()(0) << "      ";
      if (is_vector and dimension > 1) {
        file << value.real()(1) << "      ";
        if (dimension > 2) {
          file << value.real()(2) << "      ";
        }
      }
    }
    file << endl;
  }
}

void save(CMData &data, string filename) {
  save_to_file(filename, data.points, data.values, data.dimension, data.with_w,
               data.with_n, data.is_complex, data.is_vector);
}
