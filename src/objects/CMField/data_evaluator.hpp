// data_evaluator.hpp
#pragma once

#include "base_data.hpp"
//#include "cmfield.hpp"
#include "src/objects/vec.hpp"
#include "src/algorithms/spline.h"  // from tk::spline

#include <variant>
#include <stdexcept>
#include <vector>

using namespace std;

using cfloat = complex<float>;
using DataVariant = variant<
    vector<cfloat>,
    vector<vector<cfloat>>,
    vector<vector<vector<cfloat>>>,
    vector<vector<vector<vector<cfloat>>>>,
    vector<vector<vector<vector<vector<cfloat>>>>>
>;
using ResultVariant = variant<
    float,
    cfloat,
    Vec,
    complex<Vec>,
    vector<float>,                                    // For real 1D indexed arrays
    vector<cfloat>,                                   // For complex 1D indexed arrays
    vector<vector<float>>,                            // For real 2D indexed arrays (matrices)
    vector<vector<cfloat>>,                           // For complex 2D indexed arrays (matrices)
    vector<vector<vector<float>>>,                    // For real 3D indexed arrays
    vector<vector<vector<cfloat>>>,                   // For complex 3D indexed arrays
    vector<vector<vector<vector<float>>>>,            // For real 4D indexed arrays
    vector<vector<vector<vector<cfloat>>>>            // For complex 4D indexed arrays
>;


vector<Vec> invertMatrix2(vector<Vec> &matrix, int n);

inline vector<Vec> float_matrix_to_vec(vector<vector<float>> a) {
    int size = a.size();
    vector<Vec> b(size);
    for (int i = 0; i < size; i++) {
        Vec temp(a[i]);
        b[i] = temp;
    }
    return b;
}

struct DataEvaluator {
    // All values in one contiguous array, w-major like the HDF5 file:
    //   values[(w_idx * nk + spatial_idx) * block + element]
    // block is the number of values stored per grid point: the product of
    // inds for tensor fields, the vector length for vector fields, else 1.
    // Tensor elements are row-major over inds.
    vector<cfloat> values;
    int block = 1;

    vector<float> w_points;
    bool is_complex;
    bool is_vector;
    bool is_matrix;
    bool with_w;
    int dimension;
    vector<int> inds;  // Tensor indices: inds[i] = size of i-th dimension
    vector<int> mesh;
    vector<Vec> domain;
    vector<Vec> inv_domain;

    // Default constructor
    DataEvaluator()
        : is_complex(false), is_vector(false), is_matrix(false), with_w(false), dimension(1), inds({}) {}

    DataEvaluator(BaseData& f);

    ResultVariant operator()(float w);
    ResultVariant operator()(Vec point, float w = 0);
    ResultVariant operator()(Vec point, vector<int> indices, float w = 0);
    ResultVariant get_array(Vec point, float w = 0);
    ResultVariant get_array(float w);  // k-averaged for k-dependent fields

private:
    // Scattered k-points (as_mesh == false) cannot be interpolated on a grid
    bool scattered = false;
    int n_points = 1;  // spatial samples per frequency

    // (grid point, weight) pairs; a grid point indexes values in units of block
    using Taps = vector<pair<long, float>>;

    int n_spatial() const;
    int n_freq() const;
    Vec to_fractional(const Vec& point) const;
    Taps taps(const Vec* point, float w) const;   // nullptr point = k-average
    vector<cfloat> combine(const Taps& t, int offset, int count) const;
    ResultVariant scalar_result(const vector<cfloat>& v) const;
    ResultVariant tensor_result(const vector<cfloat>& v) const;
};

