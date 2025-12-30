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


inline int get_size(DataVariant &f) {
    return visit([](const auto& val) -> size_t {
        return val.size();  // Works for both vector<cfloat> and vector<vector<cfloat>>
    }, f);
}


inline vector<complex<Vec>> transform_data(DataVariant& f, int dim) {
    vector<complex<Vec>> result;

    visit([&](auto const& container) {
        using T = decay_t<decltype(container)>;

        if constexpr (is_same_v<T, vector<cfloat>>) {
            // 1D vector of cfloat
            result.reserve(container.size());
            for (auto const& val : container) {
                Vec v1(val.real()); v1.dimension = dim;
                Vec v2(val.imag()); v2.dimension = dim;
                result.emplace_back(complex<Vec>(v1, v2));
            }
        } else if constexpr (is_same_v<T, vector<vector<cfloat>>>) {
            // 2D vector of cfloat
            size_t total = 0;
            size_t size = 0;
            for (auto const& row : container) {
                size = row.size();
                total += size;
            }
            result.reserve(total);

            for (auto const& row : container) {
                Vec v1; v1.dimension = dim;
                Vec v2; v2.dimension = dim;
                int ind2 = 0;
                for (auto const& val : row) {
                    v1(ind2) = (val.real());
                    v2(ind2) = (val.imag());
                    ind2++;
                }
                result.emplace_back(complex<Vec>(v1, v2));
            }
        }
    }, f);

    return result;
}

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
    // For scalar/vector fields (rank = 0, inds = {})
    vector<complex<Vec>> data;

    // For indexed fields (rank > 0)
    // Storage layout (w-k ordering): indexed_data[w_idx][spatial_idx][indices...]
    // For example, rank=2: indexed_data_2d[w_idx][spatial_idx][i][j] where i ∈ [0, inds[0]), j ∈ [0, inds[1])
    //              rank=4: indexed_data_4d[w_idx][spatial_idx][i][j][k][l] where i ∈ [0, inds[0]), etc.
    vector<vector<vector<cfloat>>> indexed_data_1d;                         // rank=1: [w][spatial][i]
    vector<vector<vector<vector<cfloat>>>> indexed_data_2d;                 // rank=2: [w][spatial][i][j]
    vector<vector<vector<vector<vector<cfloat>>>>> indexed_data_3d;         // rank=3: [w][spatial][i][j][k]
    vector<vector<vector<vector<vector<vector<cfloat>>>>>> indexed_data_4d; // rank=4: [w][spatial][i][j][k][l]

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

    ResultVariant convert(complex<Vec>& answer);
    ResultVariant operator()(float w);
    ResultVariant operator()(Vec point, float w = 0);
    ResultVariant operator()(Vec point, vector<int> indices, float w = 0);
    ResultVariant get_array(Vec point, float w = 0);

private:
    void load_indexed_data(BaseData& f);
};

