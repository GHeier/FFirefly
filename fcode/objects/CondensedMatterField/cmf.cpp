/*
 * This file contains the implementation of the class CondensedMatterField
 * They interpolate BZ(n+1 dimension) fields
 * They have been implemented to allow for input via a file or a list of points and values.
 * They are constructed to exist over a mesh, so it is ideal for Brillouin Zone calculations.
 *
 * Author: Griffin Heier
 */
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <variant>
#include <algorithm>
#include <complex>

#include "../vec.hpp"
#include "../../algorithms/interpolate.hpp"
#include "cmf.hpp"

using namespace std;

CMF::CMF() {
    points = vector<Vec>();
    w_points = vector<float>();
    values = vector<complex<Vec>>();
    domain = vector<Vec>();
    inv_domain = vector<Vec>();
    first = Vec();
    nx = 0, ny = 0, nz = 0, nw = 0;
    wmax = 0, wmin = 0;
    dimension = 0;
    is_complex = false;
    is_vector = false;
    with_w = false;

    filled = false;
}

// Function to invert a matrix represented as vector<Vec>
vector<Vec> invertMatrix(vector<Vec>& matrix, int n) {
    // Create augmented matrix [A|I]
    vector<vector<float>> augmented(n, std::vector<float>(2 * n, 0.0f));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            augmented[i][j] = matrix[i](j);
            if (fabs(augmented[i][j]) < 1e-5) augmented[i][j] = 0.0f;
        }
        augmented[i][n + i] = 1.0f; // Identity matrix
    }

    // Gaussian Elimination
    for (size_t i = 0; i < n; ++i) {
        // Partial Pivoting
        size_t maxRow = i;
        for (size_t k = i + 1; k < n; ++k) {
            if (std::fabs(augmented[k][i]) > std::fabs(augmented[maxRow][i])) {
                maxRow = k;
            }
        }
        std::swap(augmented[i], augmented[maxRow]);

        // Check for singular matrix
        if (std::fabs(augmented[i][i]) < 1e-6) {
            throw std::runtime_error("vector<vector<float>> is singular and cannot be inverted.");
        }

        // Normalize the pivot row
        float pivot = augmented[i][i];
        for (size_t j = 0; j < 2 * n; ++j) {
            augmented[i][j] /= pivot;
        }

        // Eliminate other rows
        for (size_t k = 0; k < n; ++k) {
            if (k != i) {
                float factor = augmented[k][i];
                for (size_t j = 0; j < 2 * n; ++j) {
                    augmented[k][j] -= factor * augmented[i][j];
                }
            }
        }
    }

    // Extract the inverted matrix
    vector<Vec> inv(n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            inv[i](j) = augmented[i][n + j];
        }
    }

    return inv;
}

vector<float> matrix_multiplication(vector<Vec>& matrix, Vec& vec, int n) {
    vector<float> result(n, 0.0f);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            result[i] += matrix[i](j) * vec(j);
        }
    }
    return result;
}

Vec vec_matrix_multiplication(vector<Vec>& matrix, Vec& vec, int n) {
    Vec result; result.dimension = n;
    for (int i = 0; i < n; i++) {
        result(i) = 0;
        for (int j = 0; j < n; j++) {
            result(i) += matrix[i](j) * vec(j);
        }
    }
    return result;
}

Vec cross_product(Vec a, Vec b) {
    if (a.dimension == 3) {
        Vec result;
        result(0) = a(1) * b(2) - a(2) * b(1);
        result(1) = a(2) * b(0) - a(0) * b(2);
        result(2) = a(0) * b(1) - a(1) * b(0);
        return result;
    }
    if (a.dimension == 2) {
        Vec result;
        result(0) = a(0) * b(1) - a(1) * b(0);
        return result;
    }
    if (a.dimension == 1) {
        return a / a.norm();
    }
    throw runtime_error("Cross product only defined for 2D and 3D vectors.");
}

void CMF::get_values_for_interpolation(vector<Vec> &points, vector<float> &w_points) {
    int w_size = w_points.size();
    int section = 0;
    vector<int> section_sizes = {1, 1, 1};
    vector<Vec> domain;
    Vec first = points[0];
    this->first = first;
    int jump = 1;
    if (w_size > 1) jump = w_size;
    Vec lattice_vec = (points[jump] - first).round();

    if (dimension == 1 or (dimension == 0 and with_w)) {
        section_sizes[0] = points.size() / w_points.size();
        domain.push_back((points[points.size()-1] - first).round());
        this->domain = domain;
        inv_domain = invertMatrix(domain, dimension);
        nx = section_sizes[0];
        if (with_w) {
            wmin = *min_element(w_points.begin(), w_points.end());
            wmax = *max_element(w_points.begin(), w_points.end());
        }
        return;
    }

    for (int i = 1; i < points.size(); i += jump) {
        Vec current_vec = ((points[i] - first) / section_sizes[section]).round();
        if (cross_product(lattice_vec, current_vec).norm() < 1e-6) {
            section_sizes[section]++;
        }
        else {
            Vec p = (points[i-jump] - first).round();
            p.dimension = dimension;
            domain.push_back(p);
            jump = i;
            section++;
            section_sizes[section]++;
            if (section >= dimension) break;
            lattice_vec = (points[i] - first).round();
        }
    }
    if (w_size > 0 and domain.size() != 0) jump -= w_size;
    domain.push_back((points[points.size()-jump] - first).round());
    reverse(domain.begin(), domain.end());
    this->domain = domain;
    inv_domain = invertMatrix(domain, dimension);
    nx = section_sizes[0];
    if (with_w) nx--;
    if (dimension > 1) ny = section_sizes[1];
    if (dimension > 2) nz = section_sizes[2];
    if (dimension > 3) nw = section_sizes[3];

    if (with_w) {
        wmin = *min_element(w_points.begin(), w_points.end());
        wmax = *max_element(w_points.begin(), w_points.end());
    }
}

CMF::CMF(vector<Vec> points, vector<complex<Vec>> values, int dimension, bool with_w, bool is_complex, bool is_vector) {
    this->with_w = with_w;
    if (with_w) {
        int dim = points[0].dimension;
        dimension--;
        for (int i = 0; i < points.size(); i++) {
            if (find(w_points.begin(), w_points.end(), points[i](dim-1)) == w_points.end())
                this->w_points.push_back(points[i](dim-1));
            points[i](dim-1) = 0;
            points[i].dimension--;
        }
    }
    this->points = points;
    this->values = values;
    this->dimension = dimension;
    this->is_complex = is_complex;
    this->is_vector = is_vector;
    get_values_for_interpolation(points, w_points);
    filled = true;
}

complex<Vec> CMF::operator() (float w) {
    if (!with_w) throw runtime_error("This field does not have a w dimension.");
    return CMF_search_1d(w, w_points, values);
}

complex<Vec> CMF::operator() (Vec point, float w) {
    if (w != 0 && !with_w) throw runtime_error("This field does not have a w dimension.");
    Vec shifted = point - first;
    Vec p = vec_matrix_multiplication(inv_domain, shifted, dimension);
    p.w = w;
    if (!with_w) {
        if (dimension == 1) 
            return interpolate_1D(p.x, 0, 1, values);
        if (dimension == 2) 
            return interpolate_2D(p.x, p.y, 0, 1, 0, 1, nx, ny, values);
        if (dimension == 3) 
            return interpolate_3D(p.x, p.y, p.z, 0, 1, 0, 1, 0, 1, nx, ny, nz, values);
    }
    else {
        if (dimension == 1) 
            return CMF_search_2d(p.x, w, nx, w_points, values);
        if (dimension == 2) 
            return CMF_search_3d(p.x, p.y, w, nx, ny, w_points, values);
        if (dimension == 3) 
            return CMF_search_4d(p.x, p.y, p.z, w, nx, ny, nz, w_points, values);
    }
    throw runtime_error("Invalid dimension.");
}

complex<Vec> CMF_search_1d(float w_val, vector<float> &w_points, vector<complex<Vec>> &f) {
    // Interpolates a 1D function f(w) given on a grid between w_min and w_max
    // at the point w_val using linear interpolation
    // w_val: value at which to interpolate
    // w_min, w_max: minimum and maximum values of w
    // f: vector of function values at the grid points
    // returns: interpolated value of f(w_val)
    float w_min = w_points[0], w_max = w_points[w_points.size()-1];
    w_val = sanitize_within_bounds(w_val, w_min, w_max);
    if (w_val < w_min || w_val > w_max) throw out_of_range("w_val out of bounds");

    int i = binary_search(w_val, w_points);
    if (i == w_points.size() - 1) i--;

    float w_rel = (w_val - w_points[i]) / (w_points[i+1] - w_points[i]);
    if (w_rel < 0 || w_rel > 1) throw out_of_range("w_rel out of bounds");

    complex<Vec> result = (1 - w_rel) * f[i] + w_rel * f[i + 1];
    return result;
}

complex<Vec> CMF_search_2d(float x_val, float w_val, int nx, vector<float> &w_points, vector<complex<Vec>> &f) {
    // Interpolates a 2D function f(x, y) given on a grid between x_min and x_max
    // and y_min and y_max at the point (x_val, y_val) using bilinear interpolation
    // x_val, y_val: values at which to interpolate
    // x_min, x_max: minimum and maximum values of x
    // y_min, y_max: minimum and maximum values of y
    // f: vector of function values at the grid points
    // returns: interpolated value of f(x_val)
    float x_min = 0, x_max = 1;
    float w_min = w_points[0], w_max = w_points[w_points.size()-1];

    x_val = sanitize_within_bounds(x_val, x_min, x_max);
    w_val = sanitize_within_bounds(w_val, w_min, w_max);
    if (x_val < x_min || x_val > x_max) throw out_of_range("x_val out of bounds");
    if (w_val < w_min || w_val > w_max) throw out_of_range("w_val out of bounds");

    float dx = (x_max - x_min) / (nx - 1);

    int i = (x_val - x_min) / dx;
    int j = binary_search(w_val, w_points);
    if (i < 0 || i >= nx) throw out_of_range("i out of bounds");
    if (j < 0 || j >= w_points.size()) throw out_of_range("j out of bounds");
    if (i == nx - 1) i--;
    if (j == w_points.size() - 1) j--;

    float x_rel = (x_val - x_min) / dx - i;
    float w_rel = (w_val - w_points[j]) / (w_points[j+1] - w_points[j]);
    if (x_rel < 0 || x_rel > 1) throw out_of_range("x_rel out of bounds");
    if (w_rel < 0 || w_rel > 1) throw out_of_range("w_rel out of bounds");

    complex<Vec> result = 
        (1 - x_rel) * (1 - w_rel) * f[i*nx+j] + 
        x_rel * (1 - w_rel) * f[(i + 1)*nx+j] + 
        (1 - x_rel) * w_rel * f[i*nx+j + 1] +
        x_rel * w_rel * f[(i + 1)*nx+j + 1];

    return result;
}

complex<Vec> CMF_search_3d(float x_val, float y_val, float w_val, int nx, int ny, vector<float> &w_points, vector<complex<Vec>> &f) {
    // Interpolates a 3D function f(x, y, w) given on a grid between x_min and x_max
    // and y_min and y_max and w_min and w_max at the point (x_val, y_val, w_val) using trilinear interpolation
    // x_val, y_val, w_val: values at which to interpolate
    // x_min, x_max: minimum and maximum values of x
    // y_min, y_max: minimum and maximum values of y
    // w_min, w_max: minimum and maximum values of w
    // f: vector of function values at the grid points
    // returns: interpolated value of f(x_val, y_val, w_val)
    float x_min = 0, x_max = 1;
    float y_min = 0, y_max = 1;
    float w_min = w_points[0], w_max = w_points[w_points.size()-1];

    x_val = sanitize_within_bounds(x_val, x_min, x_max);
    y_val = sanitize_within_bounds(y_val, y_min, y_max);
    w_val = sanitize_within_bounds(w_val, w_min, w_max);
    if (x_val < x_min || x_val > x_max) throw out_of_range("x_val out of bounds");
    if (y_val < y_min || y_val > y_max) throw out_of_range("y_val out of bounds");
    if (w_val < w_min || w_val > w_max) throw out_of_range("w_val out of bounds");

    float dx = (x_max - x_min) / (nx - 1);
    float dy = (y_max - y_min) / (ny - 1);

    int i = (x_val - x_min) / dx;
    int j = (y_val - y_min) / dy;
    int k = binary_search(w_val, w_points);
    if (i < 0 || i >= nx) throw out_of_range("i out of bounds");
    if (j < 0 || j >= ny) throw out_of_range("j out of bounds");
    if (i == nx - 1) i--;
    if (j == ny - 1) j--;
    if (k == w_points.size() - 1) k--;

    float x_rel = (x_val - x_min) / dx - i;
    float y_rel = (y_val - y_min) / dy - j;
    float w_rel = (w_val - w_points[k]) / (w_points[k+1] - w_points[k]);
    if (x_rel < 0 || x_rel > 1) throw out_of_range("x_rel out of bounds");
    if (y_rel < 0 || y_rel > 1) throw out_of_range("y_rel out of bounds");
    if (w_rel < 0 || w_rel > 1) throw out_of_range("w_rel out of bounds");

    complex<Vec> result = 
        (1 - x_rel) * (1 - y_rel) * (1 - w_rel) * f[i*nx*nx+j*ny+k] + 
        x_rel * (1 - y_rel) * (1 - w_rel) * f[(i + 1)*nx*nx+j*ny+k] + 
        (1 - x_rel) * y_rel * (1 - w_rel) * f[i*nx*nx + (j + 1)*ny+k] +
        x_rel * y_rel * (1 - w_rel) * f[(i + 1)*nx*nx+(j + 1)*ny+k] +
        (1 - x_rel) * (1 - y_rel) * w_rel * f[i*nx*nx+j*ny+k + 1] +
        x_rel * (1 - y_rel) * w_rel * f[(i + 1)*nx*nx+j*ny+k + 1] +
        (1 - x_rel) * y_rel * w_rel * f[i*nx*nx + (j + 1)*ny+k + 1] +
        x_rel * y_rel * w_rel * f[(i + 1)*nx*nx+(j + 1)*ny+k + 1];

    return result;
}

complex<Vec> CMF_search_4d(float x_val, float y_val, float z_val, float w_val, 
        int nx, int ny, int nz, vector<float> &w_points, vector<complex<Vec>> &f) {
    // Interpolates a 4D function f(x, y, z, w) given on a grid between x_min and x_max
    // and y_min and y_max and z_min and z_max and w_min and w_max at the point (x_val, y_val, z_val, w_val) using trilinear interpolation
    // x_val, y_val, z_val, w_val: values at which to interpolate
    // x_min, x_max: minimum and maximum values of x
    // y_min, y_max: minimum and maximum values of y
    // z_min, z_max: minimum and maximum values of z
    // w_min, w_max: minimum and maximum values of w
    // f: vector of function values at the grid points
    // returns: interpolated value of f(x_val, y_val, z_val, w_val)
    float x_min = 0, x_max = 1;
    float y_min = 0, y_max = 1;
    float z_min = 0, z_max = 1;
    float w_min = w_points[0], w_max = w_points[w_points.size()-1];

    x_val = sanitize_within_bounds(x_val, x_min, x_max);
    y_val = sanitize_within_bounds(y_val, y_min, y_max);
    z_val = sanitize_within_bounds(z_val, z_min, z_max);
    w_val = sanitize_within_bounds(w_val, w_min, w_max);
    if (x_val < x_min || x_val > x_max) throw out_of_range("x_val out of bounds");
    if (y_val < y_min || y_val > y_max) throw out_of_range("y_val out of bounds");
    if (z_val < z_min || z_val > z_max) throw out_of_range("z_val out of bounds");
    if (w_val < w_min || w_val > w_max) throw out_of_range("w_val out of bounds");

    float dx = (x_max - x_min) / (nx - 1);
    float dy = (y_max - y_min) / (ny - 1);
    float dz = (z_max - z_min) / (nz - 1);

    int i = (x_val - x_min) / dx;
    int j = (y_val - y_min) / dy;
    int k = (z_val - z_min) / dz;
    int l = binary_search(w_val, w_points);
    if (i < 0 || i > nx) throw out_of_range("i out of bounds");
    if (j < 0 || j > ny) throw out_of_range("j out of bounds");
    if (k < 0 || k > nz) throw out_of_range("k out of bounds");
    if (i == nx - 1) i--;
    if (j == ny - 1) j--;
    if (k == nz - 1) k--;
    if (l == w_points.size() - 1) l--;

    float x_rel = (x_val - x_min) / dx - i;
    float y_rel = (y_val - y_min) / dy - j;
    float z_rel = (z_val - z_min) / dz - k;
    float w_rel = (w_val - w_points[l]) / (w_points[l+1] - w_points[l]);
    if (x_rel < 0 || x_rel > 1) throw out_of_range("x_rel out of bounds");
    if (y_rel < 0 || y_rel > 1) throw out_of_range("y_rel out of bounds");
    if (z_rel < 0 || z_rel > 1) throw out_of_range("z_rel out of bounds");
    if (w_rel < 0 || w_rel > 1) throw out_of_range("w_rel out of bounds");


    complex<Vec> result = 
        (1 - x_rel) * (1 - y_rel) * (1 - z_rel) * (1 - w_rel) * f[i*nx*nx*nx+j*ny*ny+k*nz+l] + 
        x_rel * (1 - y_rel) * (1 - z_rel) * (1 - w_rel) * f[(i + 1)*nx*nx*nx+j*ny*ny+k*nz+l] + 
        (1 - x_rel) * y_rel * (1 - z_rel) * (1 - w_rel) * f[i*nx*nx*nx + (j + 1)*ny*ny+k*nz+l] +
        x_rel * y_rel * (1 - z_rel) * (1 - w_rel) * f[(i + 1)*nx*nx*nx+(j + 1)*ny*ny+k*nz+l] +
        (1 - x_rel) * (1 - y_rel) * z_rel * (1 - w_rel) * f[i*nx*nx*nx+j*ny*ny+(k+1)*nz+l] +
        x_rel * (1 - y_rel) * z_rel * (1 - w_rel) * f[(i + 1)*nx*nx*nx+j*ny*ny+(k+1)*nz+l] +
        (1 - x_rel) * y_rel * z_rel * (1 - w_rel) * f[i*nx*nx*nx + (j + 1)*ny*ny+(k+1)*nz+l] +
        x_rel * y_rel * z_rel * (1 - w_rel) * f[(i + 1)*nx*nx*nx+(j + 1)*ny*ny+(k+1)*nz+l] +
        (1 - x_rel) * (1 - y_rel) * (1 - z_rel) * w_rel * f[i*nx*nx*nx+j*ny*ny+k*nz+l + 1] +
        x_rel * (1 - y_rel) * (1 - z_rel) * w_rel * f[(i + 1)*nx*nx*nx+j*ny*ny+k*nz+l + 1] +
        (1 - x_rel) * y_rel * (1 - z_rel) * w_rel * f[i*nx*nx*nx + (j + 1)*ny*ny+k*nz+l + 1] +
        x_rel * y_rel * (1 - z_rel) * w_rel * f[(i + 1)*nx*nx*nx+(j + 1)*ny*ny+k*nz+l + 1] +
        (1 - x_rel) * (1 - y_rel) * z_rel * w_rel * f[i*nx*nx*nx+j*ny*ny+(k+1)*nz+l + 1] +
        x_rel * (1 - y_rel) * z_rel * w_rel * f[(i + 1)*nx*nx*nx+j*ny*ny+(k+1)*nz+l + 1] +
        (1 - x_rel) * y_rel * z_rel * w_rel * f[i*nx*nx*nx + (j + 1)*ny*ny+(k+1)*nz+l + 1] +
        x_rel * y_rel * z_rel * w_rel * f[(i + 1)*nx*nx*nx+(j + 1)*ny*ny+(k+1)*nz+l + 1];

    return result;
}

CMF load_CMF_from_file(string filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        throw runtime_error("Failed to open the file.");
    }
    string line;
    int ind = 0;
    int dimension = 0;
    bool with_w = false;
    bool is_complex = false;
    bool is_vector = false;
    vector<Vec> points;
    vector<complex<Vec>> values;
    while (getline(file, line)) {
        istringstream iss(line);
        if (ind == 0) {
            string word;
            while (iss >> word) {
                if (word == "x") dimension++;
                if (word == "y") dimension++;
                if (word == "z") dimension++;
                if (word == "w") with_w = true;
                if (word == "Re(f)") is_complex = true;
                if (word == "Im(f)") is_complex = true;
                if (word == "fx") is_vector = true;
                if (word == "fy") is_vector = true;
                if (word == "fz") is_vector = true;
            }
            dimension += with_w;
            ind++;
            continue;
        }
        Vec point;
        float value;
        float ivalue;
        for (int i = 0; i < dimension; i++) {
            float temp;
            iss >> temp;
            point(i) = temp;
        }
        iss >> value;
        point.dimension = dimension;
        points.push_back(point);
        complex<Vec> vval = complex<Vec>(value, 0);
        if (is_complex) {
            iss >> ivalue;
            vval = complex<Vec>(value, ivalue);
        }
        values.push_back(vval);
    }
    return CMF(points, values, dimension, with_w, is_complex, is_vector);
}

string get_header(int dimension, bool with_w, bool is_complex, bool is_vector) {
   // Determine fheader
    string fheader = "    f    ";
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
    string header = "    x         ";
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
    if (with_w) {
        header += "   w        ";
    }

    header += fheader;
    return header;
}


void save_to_file(string filename, vector<Vec> &points, vector<complex<Vec>> &values, int dimension, bool with_w, bool is_complex, bool is_vector) {
    ofstream file(filename);
    dimension -= with_w;
    string header = get_header(dimension, with_w, is_complex, is_vector);
    file << fixed << setprecision(6);
    file << header << endl;
    for (int i = 0; i < points.size(); i++) {
        Vec point = points[i];
        complex<Vec> value = values[i];
        file << point(0) << "      ";
        if (dimension > 1) file << point(1) << "      ";
        if (dimension > 2) file << point(2) << "      ";
        if (with_w) file << point(dimension) << "      ";
        if (is_complex) {
            file << value.real()(0) << "      " << value.imag()(0) << "      ";
            if (is_vector) {
                file << value.real()(1) << "      " << value.imag()(1) << "      ";
                if (dimension > 2) {
                    file << value.real()(2) << "      " << value.imag()(2) << "      ";
                }
            }
        } else {
            file << value.real()(0) << "      ";
            if (is_vector) {
                file << value.real()(1) << "      ";
                if (dimension > 2) {
                    file << value.real()(2) << "      ";
                }
            }
        }
        file << endl;
    }
}

void combine_points_and_w(vector<Vec> &points, vector<float> &w_points) {
    for (int i = 0; i < points.size(); i++) {
        points[i].dimension++;
        points[i](points[i].dimension-1) = w_points[i % w_points.size()];
    }
}

void save_CMF_to_file(string filename, CMF &field) {
    if (field.with_w) combine_points_and_w(field.points, field.w_points);
    save_to_file(filename, field.points, field.values, field.dimension + field.with_w, field.with_w, field.is_complex, field.is_vector);
}
