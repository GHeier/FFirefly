#include "data_evaluator.hpp"
#include "cmfield.hpp"
#include "../../algorithms/interpolate.hpp"
#include <vector>
#include <complex>
#include <stdexcept>
#include <cmath>

using namespace std;
using ResultVariant = variant<
    float,
    cfloat,
    Vec,
    complex<Vec>,
    vector<float>,
    vector<cfloat>,
    vector<vector<float>>,
    vector<vector<cfloat>>
>;
using DataVariant = variant<
    vector<cfloat>,
    vector<vector<cfloat>>,
    vector<vector<vector<cfloat>>>,
    vector<vector<vector<vector<cfloat>>>>
>;

DataEvaluator::DataEvaluator(BaseData& f) {
    dimension = f.dimension;
    mesh = f.mesh;
    is_complex = f.is_complex;
    is_vector = f.is_vector;
    w_points = f.w_points;
    with_w = w_points.size() > 0;
    n_indices = f.n_indices;
    dim_indices = f.dim_indices;

    if (!f.domain.empty()) {
        domain = float_matrix_to_vec(f.domain);
        inv_domain = invertMatrix2(domain, f.dimension);
    } else {
        // Empty domain - create identity matrices
        domain.resize(dimension);
        inv_domain.resize(dimension);
        for (int i = 0; i < dimension; i++) {
            domain[i].dimension = dimension;
            inv_domain[i].dimension = dimension;
            for (int j = 0; j < dimension; j++) {
                domain[i](j) = (i == j) ? 1.0f : 0.0f;
                inv_domain[i](j) = (i == j) ? 1.0f : 0.0f;
            }
        }
    }

    if (n_indices == 0) {
        // Regular scalar/vector field
        data = transform_data(f.data, f.dimension);
    } else {
        // Indexed field - load into indexed_data structures
        load_indexed_data(f);
    }
}

void DataEvaluator::load_indexed_data(BaseData& f) {
    // Data layout expected in BaseData:
    // For n_indices=1: vector<vector<cfloat>> where outer is spatial/w, inner is index dim
    // For n_indices=2: vector<vector<vector<cfloat>>> where [spatial/w][i][j]

    if (n_indices == 1) {
        // Extract 1D indexed data
        if (auto* vec_data = std::get_if<vector<vector<cfloat>>>(&f.data)) {
            // Determine spatial dimensions
            int n_spatial = 1;
            if (f.with_k) {
                for (int m : f.mesh) n_spatial *= m;
            }
            int n_w = f.with_w ? f.w_points.size() : 1;

            // Reshape: indexed_data_1d[spatial][w][index]
            indexed_data_1d.resize(n_spatial);
            for (int s = 0; s < n_spatial; s++) {
                indexed_data_1d[s].resize(n_w);
                for (int w = 0; w < n_w; w++) {
                    int flat_idx = s * n_w + w;
                    if (flat_idx < vec_data->size()) {
                        indexed_data_1d[s][w] = (*vec_data)[flat_idx];
                    }
                }
            }
        }
    } else if (n_indices == 2) {
        // Extract 2D indexed data (matrix)
        if (auto* mat_data = std::get_if<vector<vector<vector<cfloat>>>>(&f.data)) {
            // Determine spatial dimensions
            int n_spatial = 1;
            if (f.with_k) {
                for (int m : f.mesh) n_spatial *= m;
            }
            int n_w = f.with_w ? f.w_points.size() : 1;

            // Reshape: indexed_data_2d[spatial][w][i][j]
            indexed_data_2d.resize(n_spatial);
            for (int s = 0; s < n_spatial; s++) {
                indexed_data_2d[s].resize(n_w);
                for (int w = 0; w < n_w; w++) {
                    int flat_idx = s * n_w + w;
                    if (flat_idx < mat_data->size()) {
                        indexed_data_2d[s][w] = (*mat_data)[flat_idx];
                    }
                }
            }
        }
    }
}

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

ResultVariant DataEvaluator::get_array(Vec point, float w) {
    if (n_indices == 0) {
        // No indices - just return the scalar value
        return (*this)(point, w);
    }

    // Check if data is empty
    if ((n_indices == 1 && indexed_data_1d.empty()) ||
        (n_indices == 2 && indexed_data_2d.empty())) {
        // Return empty array
        if (n_indices == 1) {
            if (is_complex) return vector<cfloat>();
            else return vector<float>();
        } else {
            if (is_complex) return vector<vector<cfloat>>();
            else return vector<vector<float>>();
        }
    }

    // Transform point to normalized coordinates
    Vec p = vec_matrix_multiplication2(inv_domain, point, dimension);
    fold_to_first_BZ2(p);

    if (n_indices == 1) {
        // Return 1D array (vector)
        vector<cfloat> result;

        if (!with_w) {
            // Spatial interpolation only
            if (dimension == 1) {
                result = interpolate_1D_vec(p.x, 0, 1, mesh[0], indexed_data_1d[0]);
            } else if (dimension == 2) {
                // Flatten indexed_data_1d for 2D interpolation
                vector<vector<cfloat>> flat_data;
                for (const auto& spatial_pt : indexed_data_1d) {
                    flat_data.push_back(spatial_pt[0]);
                }
                result = interpolate_2D_vec(p.x, p.y, 0, 1, 0, 1, mesh[0], mesh[1], flat_data);
            } else if (dimension == 3) {
                vector<vector<cfloat>> flat_data;
                for (const auto& spatial_pt : indexed_data_1d) {
                    flat_data.push_back(spatial_pt[0]);
                }
                result = interpolate_3D_vec(p.x, p.y, p.z, 0, 1, 0, 1, 0, 1, mesh[0], mesh[1], mesh[2], flat_data);
            }
        } else {
            // With frequency dimension - need to interpolate over w first, then spatial
            // For now, just return from nearest w point
            int w_idx = 0;
            for (size_t i = 0; i < w_points.size(); i++) {
                if (w >= w_points[i]) w_idx = i;
            }

            if (dimension == 1) {
                result = interpolate_1D_vec(p.x, 0, 1, mesh[0], indexed_data_1d[0]);
            } else if (dimension == 2) {
                vector<vector<cfloat>> flat_data;
                for (const auto& spatial_pt : indexed_data_1d) {
                    flat_data.push_back(spatial_pt[w_idx]);
                }
                result = interpolate_2D_vec(p.x, p.y, 0, 1, 0, 1, mesh[0], mesh[1], flat_data);
            } else if (dimension == 3) {
                vector<vector<cfloat>> flat_data;
                for (const auto& spatial_pt : indexed_data_1d) {
                    flat_data.push_back(spatial_pt[w_idx]);
                }
                result = interpolate_3D_vec(p.x, p.y, p.z, 0, 1, 0, 1, 0, 1, mesh[0], mesh[1], mesh[2], flat_data);
            }
        }

        // Convert to real if needed
        if (!is_complex) {
            vector<float> real_result(result.size());
            for (size_t i = 0; i < result.size(); i++) {
                real_result[i] = result[i].real();
            }
            return real_result;
        }
        return result;

    } else if (n_indices == 2) {
        // Return 2D array (matrix)
        vector<vector<cfloat>> result;

        if (!with_w) {
            // Spatial interpolation only - extract w=0 data for all spatial points
            if (dimension == 1) {
                vector<vector<vector<cfloat>>> spatial_data;
                for (const auto& spatial_pt : indexed_data_2d) {
                    spatial_data.push_back(spatial_pt[0]);
                }
                result = interpolate_1D_mat(p.x, 0, 1, mesh[0], spatial_data);
            } else if (dimension == 2) {
                vector<vector<vector<cfloat>>> flat_data;
                for (const auto& spatial_pt : indexed_data_2d) {
                    flat_data.push_back(spatial_pt[0]);
                }
                result = interpolate_2D_mat(p.x, p.y, 0, 1, 0, 1, mesh[0], mesh[1], flat_data);
            } else if (dimension == 3) {
                vector<vector<vector<cfloat>>> flat_data;
                for (const auto& spatial_pt : indexed_data_2d) {
                    flat_data.push_back(spatial_pt[0]);
                }
                result = interpolate_3D_mat(p.x, p.y, p.z, 0, 1, 0, 1, 0, 1, mesh[0], mesh[1], mesh[2], flat_data);
            }
        } else {
            // With frequency dimension - interpolate in frequency
            // Find the two w points to interpolate between
            int w_idx = 0;
            for (size_t i = 0; i < w_points.size(); i++) {
                if (w >= w_points[i]) w_idx = i;
            }
            if (w_idx >= w_points.size() - 1) w_idx = w_points.size() - 2;

            // Linear interpolation weight in frequency
            float w_weight = 0.0;
            if (w_points.size() > 1 && w_idx < w_points.size() - 1) {
                w_weight = (w - w_points[w_idx]) / (w_points[w_idx + 1] - w_points[w_idx]);
            }

            // Get results at both w points
            vector<vector<cfloat>> result_w0, result_w1;

            if (dimension == 1) {
                vector<vector<vector<cfloat>>> spatial_data_w0, spatial_data_w1;
                for (const auto& spatial_pt : indexed_data_2d) {
                    spatial_data_w0.push_back(spatial_pt[w_idx]);
                    spatial_data_w1.push_back(spatial_pt[w_idx + 1]);
                }
                result_w0 = interpolate_1D_mat(p.x, 0, 1, mesh[0], spatial_data_w0);
                result_w1 = interpolate_1D_mat(p.x, 0, 1, mesh[0], spatial_data_w1);
            } else if (dimension == 2) {
                vector<vector<vector<cfloat>>> flat_data_w0, flat_data_w1;
                for (const auto& spatial_pt : indexed_data_2d) {
                    flat_data_w0.push_back(spatial_pt[w_idx]);
                    flat_data_w1.push_back(spatial_pt[w_idx + 1]);
                }
                result_w0 = interpolate_2D_mat(p.x, p.y, 0, 1, 0, 1, mesh[0], mesh[1], flat_data_w0);
                result_w1 = interpolate_2D_mat(p.x, p.y, 0, 1, 0, 1, mesh[0], mesh[1], flat_data_w1);
            } else if (dimension == 3) {
                vector<vector<vector<cfloat>>> flat_data_w0, flat_data_w1;
                for (const auto& spatial_pt : indexed_data_2d) {
                    flat_data_w0.push_back(spatial_pt[w_idx]);
                    flat_data_w1.push_back(spatial_pt[w_idx + 1]);
                }
                result_w0 = interpolate_3D_mat(p.x, p.y, p.z, 0, 1, 0, 1, 0, 1, mesh[0], mesh[1], mesh[2], flat_data_w0);
                result_w1 = interpolate_3D_mat(p.x, p.y, p.z, 0, 1, 0, 1, 0, 1, mesh[0], mesh[1], mesh[2], flat_data_w1);
            }

            // Interpolate between the two frequency points
            int rows = result_w0.size();
            int cols = result_w0.empty() ? 0 : result_w0[0].size();
            result.resize(rows, vector<cfloat>(cols));
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    result[i][j] = (1.0f - w_weight) * result_w0[i][j] + w_weight * result_w1[i][j];
                }
            }
        }

        // Convert to real if needed
        if (!is_complex) {
            if (result.empty() || result[0].empty()) {
                return vector<vector<float>>();
            }
            vector<vector<float>> real_result(result.size(), vector<float>(result[0].size()));
            for (size_t i = 0; i < result.size(); i++) {
                for (size_t j = 0; j < result[i].size(); j++) {
                    real_result[i][j] = result[i][j].real();
                }
            }
            return real_result;
        }
        return result;
    } else {
        throw std::runtime_error("n_indices > 2 not yet supported");
    }
}

ResultVariant DataEvaluator::operator()(Vec point, vector<int> indices, float w) {
    if (n_indices == 0) {
        throw std::runtime_error("Indexed operator called on non-indexed field");
    }

    if (indices.size() != n_indices) {
        throw std::runtime_error("Number of indices doesn't match n_indices");
    }

    // Calculate the linear index from multi-dimensional indices
    int linear_idx = 0;
    int multiplier = 1;
    for (int i = n_indices - 1; i >= 0; i--) {
        if (indices[i] < 0 || indices[i] >= dim_indices) {
            throw std::runtime_error("Index out of bounds");
        }
        linear_idx += indices[i] * multiplier;
        multiplier *= dim_indices;
    }

    // For now, just evaluate the point normally and return a placeholder
    // Full implementation requires storing data[k_idx][w_idx][index_0][index_1]...
    auto base_result = (*this)(point, w);

    // Convert to array result
    if (is_complex) {
        vector<cfloat> result(1);
        if (auto* val = std::get_if<cfloat>(&base_result)) {
            result[0] = *val;
        }
        return result;
    } else {
        vector<float> result(1);
        if (auto* val = std::get_if<float>(&base_result)) {
            result[0] = *val;
        }
        return result;
    }
}

