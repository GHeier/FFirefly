#include "field_evaluator.hpp"
#include "cmfield.hpp"
#include <vector>
#include <complex>
#include <stdexcept>

using namespace std;
using ResultVariant = variant<
    float,
    cfloat,
    Vec,
    complex<Vec>
>;
using DataVariant = variant<
    vector<cfloat>,         
    vector<vector<cfloat>>  
>;

ResultVariant FieldEvaluator::convert(complex<Vec>& answer) {
    if (is_complex) {
        if (is_vector) {
            return answer;
        }
        return complex<float>(answer.real()(0), answer.imag()(0));
    }
    return float(answer.real()(0));
}

ResultVariant FieldEvaluator::operator()(float w) {
    complex<Vec> answer = CMF_search_1d(w, w_points, data);
    return convert(answer);
}

ResultVariant FieldEvaluator::operator()(Vec point, float w) {
    Vec shifted = point - first;
    Vec p = vec_matrix_multiplication(inv_domain, shifted, data.dimension);
    fold_to_first_BZ(p);
    p.w = w;
    if (!data.with_w) {
        if (data.dimension == 1)
            return interpolate_1D(p.x, 0, 1, values[0]);
        if (data.dimension == 2)
            return interpolate_2D(p.x, p.y, 0, 1, 0, 1, nx, ny, values[0]);
        if (data.dimension == 3)
            return interpolate_3D(p.x, p.y, p.z, 0, 1, 0, 1, 0, 1, nx, ny, nz,
                                  values[0]);
    } else {
        if (data.dimension == 1)
            return CMF_search_2d(p.x, w, nx, data.w_points, values[0]);
        if (data.dimension == 2)
            return CMF_search_3d(p.x, p.y, w, nx, ny, data.w_points, values[0]);
        if (data.dimension == 3)
            return CMF_search_4d(p.x, p.y, p.z, w, nx, ny, nz, data.w_points,
                                 values[0]);
    }
    throw runtime_error("Invalid dimension.");
}

