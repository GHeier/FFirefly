#include "data_evaluator.hpp"
#include "cmfield.hpp"
#include "../../algorithms/interpolate.hpp"
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

ResultVariant DataEvaluator::convert(complex<Vec>& answer) {
    if (is_complex) {
        if (is_vector) {
            return answer;
        }
        return complex<float>(answer.real()(0), answer.imag()(0));
    }
    return float(answer.real()(0));
}

vector<Vec> invertMatrix2(vector<Vec> &matrix, int n) {
    // Create augmented matrix [A|I]
    vector<vector<float>> augmented(n, std::vector<float>(2 * n, 0.0f));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            augmented[i][j] = matrix[i](j);
            if (fabs(augmented[i][j]) < 1e-5)
                augmented[i][j] = 0.0f;
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
            throw std::runtime_error(
                "vector<vector<float>> is singular and cannot be inverted.");
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

double fold2(double x) {
    // return x;
    double decimal = x - std::floor(x);
    return decimal;
}

void fold_to_first_BZ2(Vec &p) {
    p.x = fold2(p.x);
    p.y = fold2(p.y);
    p.z = fold2(p.z);
}

ResultVariant DataEvaluator::operator()(float w) {
    complex<Vec> answer = CMF_search_1d(w, w_points, data);
    return convert(answer);
}

Vec vec_matrix_multiplication2(vector<Vec> &matrix, Vec &vec, int n) {
    Vec result;
    result.dimension = n;
    for (int i = 0; i < n; i++) {
        result(i) = 0;
        for (int j = 0; j < n; j++) {
            result(i) += matrix[i](j) * vec(j);
        }
    }
    return result;
}

ResultVariant DataEvaluator::operator()(Vec point, float w) {
    Vec p = vec_matrix_multiplication2(inv_domain, point, dimension);
    fold_to_first_BZ2(p);
    p.w = w;

    complex<Vec> ans;
    if (!with_w) {
        if (dimension == 1)
            ans = interpolate_1D(p.x, 0, 1, data);
        else if (dimension == 2)
            ans = interpolate_2D(p.x, p.y, 0, 1, 0, 1, mesh[0], mesh[1], data);
        else if (dimension == 3)
            ans = interpolate_3D(p.x, p.y, p.z, 0, 1, 0, 1, 0, 1, mesh[0], mesh[1], mesh[2],
                                  data);
    } else {
        if (dimension == 1)
            ans = CMF_search_2d(p.x, w, mesh[0], w_points, data);
        else if (dimension == 2)
            ans = CMF_search_3d(p.x, p.y, w, mesh[0], mesh[1], w_points, data);
        else if (dimension == 3)
            ans = CMF_search_4d(p.x, p.y, p.z, w, mesh[0], mesh[1], mesh[2], w_points,
                                 data);
    }
    return convert(ans);
}

