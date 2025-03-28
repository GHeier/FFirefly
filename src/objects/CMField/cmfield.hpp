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

#include "../vec.hpp"
#include "../CMData/cmdata.hpp"

using namespace std;
class CMF {
    public:
        CMData data;
        vector<vector<complex<Vec>>> values;

        vector<Vec> domain;
        vector<Vec> inv_domain;
        Vec first;
        float xmin, xmax, ymin, ymax, zmin, zmax, wmin, wmax;
        int nx = 0, ny = 0, nz = 0, nw = 0;

        CMF();
        CMF(CMData &data);
        CMF(vector<Vec> points, vector<complex<Vec>> values, int coords_dim, bool w_col, bool n_col, bool is_complex, bool is_vector);

        void find_domain(CMData &data);
        void find_domain_1d(CMData &data);
        complex<Vec> operator() (float w);
        complex<Vec> operator() (int n, float w);
        complex<Vec> operator() (Vec point, float w = 0);
        complex<Vec> operator() (int n, Vec point, float w = 0);
};


vector<float> matrix_multiplication(vector<Vec>& matrix, Vec& vec, int n);
Vec vec_matrix_multiplication(vector<float>& matrix, Vec& vec, int n);
bool domain_vec_found(float a[4], float dx, float dy, float dz, float dw);
vector<Vec> invertMatrix(vector<Vec>& matrix, int n);
complex<Vec> CMF_search_1d(float w_val, vector<float> &w_points, vector<complex<Vec>> &f);
complex<Vec> CMF_search_2d(float x_val, float w_val, int nx, vector<float> &w_points, vector<complex<Vec>> &f);
complex<Vec> CMF_search_3d(float x_val, float y_val, float w_val, int nx, int ny, vector<float> &w_points, vector<complex<Vec>> &f);
complex<Vec> CMF_search_4d(float x_val, float y_val, float z_val, float w_val, int nx, int ny, int nz, vector<float> &w_points, vector<complex<Vec>> &f);

void empty_CMData(CMData &data);
void make_values_2d(vector<vector<complex<Vec>>> &values, CMData &data);
CMF load_CMF(string filename);
void save_CMF(string filename, CMF &field);
