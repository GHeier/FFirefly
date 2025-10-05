// data_evaluator.hpp
#pragma once

#include "base_data.hpp"
#include "cmfield.hpp"
#include "../vec.hpp"
#include "../../algorithms/spline.h"  // from tk::spline

#include <variant>
#include <stdexcept>
#include <vector>

using namespace std;

using cfloat = complex<float>;
using DataVariant = variant<
    vector<cfloat>,
    vector<vector<cfloat>>,
    vector<vector<vector<cfloat>>>,
    vector<vector<vector<vector<cfloat>>>>
>;
using ResultVariant = variant<
    float,
    cfloat,
    Vec,
    complex<Vec>,
    vector<float>,           // For real indexed arrays
    vector<cfloat>,          // For complex indexed arrays
    vector<vector<float>>,   // For real multi-dimensional indexed arrays
    vector<vector<cfloat>>   // For complex multi-dimensional indexed arrays
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

            int ind1 = 0;
            for (auto const& row : container) {
                Vec v1; v1.dimension = dim;
                Vec v2; v2.dimension = dim;
                int ind2 = 0;
                for (auto const& val : row) {
                    v1(ind2) = (val.real());
                    v2(ind2) = (val.imag());
                    ind2++;
                }
                result[ind1] = (complex<Vec>(v1, v2));
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
    // For scalar/vector fields (n_indices = 0)
    vector<complex<Vec>> data;

    // For indexed fields (n_indices > 0)
    // Storage layout: indexed_data[spatial_idx][w_idx][flat_index]
    // where flat_index = i0 + i1*dim_indices + i2*dim_indices^2 + ...
    vector<vector<vector<cfloat>>> indexed_data_1d;  // For 1D indexed (vectors)
    vector<vector<vector<vector<cfloat>>>> indexed_data_2d;  // For 2D indexed (matrices)

    vector<float> w_points;
    bool is_complex;
    bool is_vector;
    bool is_matrix;
    bool with_w;
    int dimension;
    int n_indices;
    int dim_indices;
    vector<int> mesh;
    vector<Vec> domain;
    vector<Vec> inv_domain;

    // Default constructor
    DataEvaluator()
        : is_complex(false), is_vector(false), is_matrix(false), with_w(false), dimension(1), n_indices(0), dim_indices(1) {}

    DataEvaluator(BaseData& f);

    ResultVariant convert(complex<Vec>& answer);
    ResultVariant operator()(float w);
    ResultVariant operator()(Vec point, float w = 0);
    ResultVariant operator()(Vec point, vector<int> indices, float w = 0);
    ResultVariant get_array(Vec point, float w = 0);

private:
    void load_indexed_data(BaseData& f);
};

