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

struct FieldEvaluator {
    vector<complex<Vec>> data;
    vector<float> w_points;
    bool is_complex;
    bool is_vector;
    Vec first(-0.5, -0.5, -0.5, 0, 0, 1);

    FieldEvaluator(BaseField& f) {
        data = transform_data(f.data, f.dimension);
        is_complex = f.is_complex;
        is_vector = f.is_vector;
        w_points = f.w_points;
        first.dimension = f.dimension;
    }

    ResultVariant convert(complex<Vec>& answer);
    ResultVariant operator()(float w);
};

