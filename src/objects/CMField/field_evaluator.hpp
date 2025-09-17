// field_evaluator.hpp
#pragma once

#include "base_field.hpp"
#include "cmfield.hpp"
#include "../../algorithms/spline.h"  // from tk::spline

#include <variant>
#include <stdexcept>
#include <vector>

using namespace std;

using cfloat = complex<float>;
using DataVariant = variant<
    vector<cfloat>,         
    vector<vector<cfloat>>  
>;
using ResultVariant = variant<
    float,
    cfloat,
    Vec,
    complex<Vec>
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
                Vec v(val.real(), val.imag(), 0.0f, 0.0f, 0.0f, dim, 1);
                result.emplace_back(v);  // wraps Vec into complex<Vec>
            }
        } else if constexpr (is_same_v<T, vector<vector<cfloat>>>) {
            // 2D vector of cfloat
            size_t total = 0;
            for (auto const& row : container) total += row.size();
            result.reserve(total);

            for (auto const& row : container) {
                for (auto const& val : row) {
                    Vec v(val.real(), val.imag(), 0.0f, 0.0f, 0.0f, dim, 1);
                    result.emplace_back(v);
                }
            }
        }
    }, f);

    return result;
}

ResultVariant Field_search_1d(float w_val, const vector<float>& w_points, const vector<complex<Vec>>& f);

struct FieldEvaluator {
    vector<complex<Vec>> data;
    vector<float> w_points;
    bool is_complex;
    bool is_vector;

    FieldEvaluator(BaseField& f) {
        data = transform_data(f.data, f.dimension);
        is_complex = f.is_complex;
        is_vector = f.is_vector;
        w_points = f.w_points;
    }

    ResultVariant operator()(float w);
};

