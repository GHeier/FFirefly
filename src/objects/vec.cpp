/**
 * This file contains the implementation of the Vec class
 *
 * Author: Griffin Heier
 */

#include "vec.hpp"
#include <complex>
#include <math.h>
#include <string>
#include <vector>

Vec::Vec() {
  dimension = 3;
  area = 0;
  x = 0;
  y = 0;
  z = 0;
  w = 0;
  n = 1;
}

Vec::Vec(float _x, float _y, float _z, float _w, float _area, int _dimension,
         int _n) {
  dimension = _dimension;
  area = _area;
  x = _x;
  y = _y;
  z = _z;
  w = _w;
  n = _n;
}

Vec::Vec(const float *points, int len) {
  dimension = len;
  area = 0;
  if (dimension > 0)
    x = points[0];
  if (dimension > 1)
    y = points[1];
  if (dimension > 2)
    z = points[2];
  if (dimension > 3)
    w = points[3];
  n = 1;
}

Vec::Vec(vector<float> input) {
  dimension = input.size();
  area = 0;
  if (dimension > 0)
    x = input[0];
  if (dimension > 1)
    y = input[1];
  if (dimension > 2)
    z = input[2];
  if (dimension > 3)
    w = input[3];
  n = 1;
}

float &Vec::operator()(int i) {
  if (i == 0)
    return x;
  if (i == 1)
    return y;
  if (i == 2)
    return z;
  if (i == 3)
    return w;
  printf("Invalid index of %d for Vec\n", i);
  exit(1);
}

Vec Vec::round(int precision) {
  float r = pow(16, precision);
  Vec result(std::round(x * r) / r, std::round(y * r) / r,
             std::round(z * r) / r, std::round(w * r) / r, area, dimension, n);
  return result;
}

Vec string_to_vec(string str) {
  int start = 0;
  int end = str.find(" ");
  Vec result;
  int iter = 0;
  while (end != -1) {
    string temp = str.substr(start, end - start);
    float val = std::stod(temp);
    if (iter == 0)
      result.x = val;
    if (iter == 1)
      result.y = val;
    if (iter == 2)
      result.z = val;
    if (iter == 3)
      result.w = val;
    start = end + 1;
    end = str.find(" ", start);
    iter++;
  }
  return result;
}

vector<float> unpack_string(string str) {
  vector<string> result(4);
  int start = 0;
  int end = str.find(" ");
  int iter = 0;
  while (end != -1) {
    string temp = str.substr(start, end - start);
    result[iter] = temp;
    start = end + 1;
    end = str.find(" ", start);
    iter++;
  }
  result[3] = str.substr(start, str.length() - start);
  vector<float> result_float(4);
  for (int i = 0; i < 4; i++) {
    if (result[i] == "") {
      result_float[i] = 0;
      continue;
    }
    result_float[i] = std::stod(result[i]);
  }
  return result_float;
}

string vec_to_string(Vec k) {
  string result = to_string(k.x) + " " + to_string(k.y) + " " + to_string(k.z) +
                  " " + to_string(k.w);
  return result;
}

Vec operator+(const Vec &k, const Vec &q) {
  return Vec(k.x + q.x, k.y + q.y, k.z + q.z, k.w + q.w, k.area, k.dimension,
             k.n);
}

Vec operator-(const Vec &k, const Vec &q) {
  return Vec(k.x - q.x, k.y - q.y, k.z - q.z, k.w - q.w, k.area, k.dimension,
             k.n);
}

Vec operator*(const Vec &input, int multiple) {
  return Vec(input.x * multiple, input.y * multiple, input.z * multiple,
             input.w * multiple, input.area, input.dimension, input.n);
}

Vec operator*(int multiple, const Vec &input) {
  return Vec(input.x * multiple, input.y * multiple, input.z * multiple,
             input.w * multiple, input.area, input.dimension, input.n);
}

Vec operator*(const Vec &input, double multiple) {
  return Vec(input.x * multiple, input.y * multiple, input.z * multiple,
             input.w * multiple, input.area, input.dimension, input.n);
}

Vec operator*(double multiple, const Vec &input) {
  return Vec(input.x * multiple, input.y * multiple, input.z * multiple,
             input.w * multiple, input.area, input.dimension, input.n);
}

Vec operator*(vector<vector<double>> &left, Vec right) {
  int dim = left.size() > right.dimension ? right.dimension : left.size();
  Vec result;
  result.dimension = right.dimension;
  result.n = right.n;
  result.area = right.area;
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      result(i) += left[i][j] * right(j);
    }
  }
  return result;
}

Vec operator/(const Vec &input, double multiple) {
  return Vec(input.x / multiple, input.y / multiple, input.z / multiple,
             input.w / multiple, input.area, input.dimension, input.n);
}

Vec operator/(double multiple, const Vec &input) {
  return Vec(input.x / multiple, input.y / multiple, input.z / multiple,
             input.w / multiple, input.area, input.dimension, input.n);
}

Vec operator*(const Vec &input, float multiple) {
  return Vec(input.x * multiple, input.y * multiple, input.z * multiple,
             input.w * multiple, input.area, input.dimension, input.n);
}

Vec operator*(float multiple, const Vec &input) {
  return Vec(input.x * multiple, input.y * multiple, input.z * multiple,
             input.w * multiple, input.area, input.dimension, input.n);
}

float operator*(const Vec &left, const Vec &right) {
  return left.x * right.x + left.y * right.y + left.z * right.z +
         left.w * right.w;
}

Vec operator*(vector<vector<float>> &left, Vec right) {
  int dim = left.size() > right.dimension ? right.dimension : left.size();
  Vec result;
  result.dimension = right.dimension;
  result.n = right.n;
  result.area = right.area;
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      result(i) += left[i][j] * right(j);
    }
  }
  return result;
}

Vec operator/(const Vec &input, float multiple) {
  return Vec(input.x / multiple, input.y / multiple, input.z / multiple,
             input.w / multiple, input.area, input.dimension, input.n);
}

Vec operator/(float multiple, const Vec &input) {
  return Vec(input.x / multiple, input.y / multiple, input.z / multiple,
             input.w / multiple, input.area, input.dimension, input.n);
}

Vec operator/(const Vec &input, int multiple) {
  return Vec(input.x / multiple, input.y / multiple, input.z / multiple,
             input.w / multiple, input.area, input.dimension, input.n);
}

Vec operator/(int multiple, const Vec &input) {
  return Vec(input.x / multiple, input.y / multiple, input.z / multiple,
             input.w / multiple, input.area, input.dimension, input.n);
}

complex<Vec> operator*(const float &left, const complex<Vec> &right) {
  return complex<Vec>(left * right.real(), left * right.imag());
}

complex<Vec> operator*(const complex<Vec> &left, const float &right) {
  return complex<Vec>(left.real() * right, left.imag() * right);
}

complex<Vec> operator/(const complex<Vec> &left, const float &right) {
  return complex<Vec>(left.real() / right, left.imag() / right);
}

complex<Vec> operator/(const float &left, const complex<Vec> &right) {
  return complex<Vec>(left / right.real(), left / right.imag());
}

complex<Vec> operator+(const complex<Vec> &left, const complex<Vec> &right) {
  return complex<Vec>(left.real() + right.real(), left.imag() + right.imag());
}

complex<Vec> operator+(const complex<Vec> &left, const Vec &right) {
  return complex<Vec>(left.real() + right, left.imag());
}

complex<Vec> operator+(const Vec &left, const complex<Vec> &right) {
  return complex<Vec>(left + right.real(), right.imag());
}

complex<Vec> operator-(const complex<Vec> &left, const complex<Vec> &right) {
  return complex<Vec>(left.real() - right.real(), left.imag() - right.imag());
}

complex<Vec> operator-(const complex<Vec> &left, const Vec &right) {
  return complex<Vec>(left.real() - right, left.imag());
}

complex<Vec> operator-(const Vec &left, const complex<Vec> &right) {
  return complex<Vec>(left - right.real(), float(-1.0) * right.imag());
}

float Vec::norm() {
  if (dimension == 1)
    return sqrt(x * x);
  if (dimension == 2)
    return sqrt(x * x + y * y);
  if (dimension == 3)
    return sqrt(x * x + y * y + z * z);
  return sqrt(x * x + y * y + z * z + w * w);
}

bool operator==(const Vec &k, const Vec &q) {
  return k.x == q.x && k.y == q.y && k.z == q.z && k.w == q.w;
}

bool operator<(const Vec &left, const Vec &right) { return left.w < right.w; }

std::ostream &operator<<(std::ostream &os, const Vec &k) {
  if (k.dimension == 1)
    os << k.x;
  if (k.dimension == 2)
    os << k.x << " " << k.y;
  if (k.dimension == 3)
    os << k.x << " " << k.y << " " << k.z;
  if (k.dimension == 4)
    os << k.x << " " << k.y << " " << k.z << " " << k.w;
  return os;
}
