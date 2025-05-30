/**
 * @class Vec
 * @brief A class to represent a vector in 3D space.
 *
 * This is a vector that could be on the Fermi Surface, and so makes up a small
 * triangle on that surface. As such, it has an area and a "frequency"
 * associated with it. The frequency is the energy associated with the vector,
 * as it can also be off of the Fermi Surface, and so will have a nonzero
 * energy.
 */
#pragma once

#include <complex>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

class Vec {
  public:
    float x;
    float y;
    float z;
    float w;
    float area = 0;
    int dimension = 3;
    int n = 1;

    Vec();
    explicit Vec(float _x, float _y = 0, float _z = 0, float _w = 0, float _area = 0, int _dimension = 3, int _n = 1);
    Vec(vector<float> input);
    Vec(const float *points, int len);
    float &operator()(int i);
    Vec round(int precision = 4);
    float norm();
};

/// A type alias for Vec when converting a vector to and from a string for a
/// hash map.
Vec string_to_vec(string str);
vector<float> unpack_string(string str);
/// A type alias for Vec when converting a vector to and from a string for a
/// hash map.
string vec_to_string(Vec k);
Vec operator+(const Vec &k, const Vec &q);
Vec operator-(const Vec &k, const Vec &q);
Vec operator*(int multiple, const Vec &input);
Vec operator*(const Vec &input, int multiple);
Vec operator*(float multiple, const Vec &input);
Vec operator*(const Vec &input, float multiple);
Vec operator*(double multiple, const Vec &input);
Vec operator*(const Vec &input, double multiple);
/// Dot product of two vectors.
float operator*(const Vec &left, const Vec &right);
Vec operator*(vector<vector<float>> &left, Vec right);
Vec operator/(int multiple, const Vec &input);
Vec operator/(const Vec &input, int multiple);
Vec operator/(float multiple, const Vec &input);
Vec operator/(const Vec &input, float multiple);
Vec operator/(double multiple, const Vec &input);
Vec operator/(const Vec &input, double multiple);

complex<Vec> operator*(const float &left, const complex<Vec> &right);
complex<Vec> operator*(const complex<Vec> &left, const float &right);
complex<Vec> operator/(const complex<Vec> &left, const float &right);
complex<Vec> operator/(const float &left, const complex<Vec> &right);
complex<Vec> operator+(const complex<Vec> &left, const complex<Vec> &right);
complex<Vec> operator+(const complex<Vec> &left, const Vec &right);
complex<Vec> operator+(const Vec &left, const complex<Vec> &right);
complex<Vec> operator-(const complex<Vec> &left, const complex<Vec> &right);
complex<Vec> operator-(const complex<Vec> &left, const Vec &right);
complex<Vec> operator-(const Vec &left, const complex<Vec> &right);

bool operator==(const Vec &k, const Vec &q);
/// Comparison operator for Vec, sorting them by area
bool operator<(const Vec &left, const Vec &right);
ostream &operator<<(std::ostream &os, const Vec &k);
