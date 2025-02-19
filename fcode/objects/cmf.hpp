/*
 * This file contains the definition of the ScalarField and VectorField classes.
 * These classes are used to represent scalar and vector fields, respectively.
 * They are implemented over a mesh of points and values, and can be used to interpolate
 * The mesh is interpolated automatically when the field is created. Ideal for BZ calculations
 *
 * Author: Griffin Heier
 */
#pragma once

#include <string>
#include <complex>
#include <vector>

#include "vec.hpp"

using namespace std;
class CMF {
    public:
        vector<Vec> points;
        vector<float> w_points;
        vector<complex<Vec>> values;
        vector<Vec> ivalues;
        vector<Vec> domain;
        vector<Vec> inv_domain;
        Vec first;
        float xmin, xmax, ymin, ymax, zmin, zmax, wmin, wmax;
        int nx, ny, nz, nw;
        int dimension;
        bool is_complex;
        bool is_vector;
        bool with_w;
        bool filled;

        CMF();
        CMF(vector<Vec> points, vector<complex<Vec>> values, int dimension, bool with_w, bool is_complex, bool is_vector);

        void get_values_for_interpolation(vector<Vec> &points, vector<float> &w_points);
        complex<Vec> operator() (Vec point, float w = 0);
        complex<Vec> operator() (float w);
};


vector<float> matrix_multiplication(vector<Vec>& matrix, Vec& vec, int n);
Vec vec_matrix_multiplication(vector<float>& matrix, Vec& vec, int n);
bool domain_vec_found(float a[4], float dx, float dy, float dz, float dw);
vector<Vec> invertMatrix(vector<Vec>& matrix, int n);
complex<Vec> CMF_search_1d(float w_val, vector<float> &w_points, vector<complex<Vec>> &f);
complex<Vec> CMF_search_2d(float x_val, float w_val, int nx, vector<float> &w_points, vector<complex<Vec>> &f);
complex<Vec> CMF_search_3d(float x_val, float y_val, float w_val, int nx, int ny, vector<float> &w_points, vector<complex<Vec>> &f);
complex<Vec> CMF_search_4d(float x_val, float y_val, float z_val, float w_val, int nx, int ny, int nz, vector<float> &w_points, vector<complex<Vec>> &f);

CMF load_CMF_from_file(string filename);
void save_to_file(string filename, vector<Vec> &points, vector<complex<Vec>> &values, int dimension, bool with_w, bool is_complex, bool is_vector);
void save_CMF_to_file(string filename, CMF &field);
