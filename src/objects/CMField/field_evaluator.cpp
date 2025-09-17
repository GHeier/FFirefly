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

ResultVariant FieldEvaluator::operator()(float w) {
    complex<Vec> answer = CMF_search_1d(w, w_points, data);
    if (is_complex) {
        if (is_vector) {
            return answer;
        }
        return complex<float>(answer.real()(0), answer.imag()(0));
    }
    return float(answer.real()(0));
}

// cmfield.hpp or a dedicated interpolation header
// Placeholder for binary search utility — assumes sorted `w_points`
int binary_search(float w, const vector<float>& w_points) {
    auto it = upper_bound(w_points.begin(), w_points.end(), w);
    int idx = max(0, int(it - w_points.begin()) - 1);
    return idx;
}

// Optional: clamp to [min, max] if needed
float sanitize_within_bounds(float w, float min_w, float max_w) {
    return clamp(w, min_w, max_w); // Requires C++17
}
