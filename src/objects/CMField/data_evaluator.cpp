#include "data_evaluator.hpp"
#include "field_funcs.hpp"
#include "src/algorithms/interpolate.hpp"
#include <vector>
#include <complex>
#include <stdexcept>
#include <cmath>

using namespace std;
// Note: ResultVariant and DataVariant are defined in data_evaluator.hpp

// Number of spatial samples stored per frequency
static int spatial_count(const BaseData& f) {
    if (!f.mesh.empty()) {
        int n = 1;
        for (int m : f.mesh) n *= m;
        return n;
    }
    if (f.with_k && !f.points.empty()) return f.points.size();
    return 1;
}

// A domain is usable if it has at least `dim` rows of at least `dim` entries.
// The legacy default {{}} (one empty row) counts as "no domain".
static bool domain_usable(const vector<vector<float>>& domain, int dim) {
    if (dim <= 0 || (int)domain.size() < dim) return false;
    for (int i = 0; i < dim; i++)
        if ((int)domain[i].size() < dim) return false;
    return true;
}

static bool domain_absent(const vector<vector<float>>& domain) {
    for (auto& row : domain)
        if (!row.empty()) return false;
    return true;
}

DataEvaluator::DataEvaluator(BaseData& f) {
    dimension = f.dimension;
    mesh = f.mesh;
    is_complex = f.is_complex;
    is_vector = f.is_vector;
    w_points = f.w_points;
    with_w = w_points.size() > 0;
    inds = f.inds;
    scattered = f.with_k && !f.as_mesh && f.mesh.empty() && !f.points.empty();

    if (domain_usable(f.domain, dimension)) {
        // Validate domain before conversion and inversion
        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                if (std::isnan(f.domain[i][j]) || std::isinf(f.domain[i][j])) {
                    throw std::runtime_error("Domain matrix contains NaN or Inf in DataEvaluator constructor at ["
                                           + std::to_string(i) + "][" + std::to_string(j) + "]. "
                                           "This indicates memory corruption before DataEvaluator initialization.");
                }
            }
        }

        domain = float_matrix_to_vec(f.domain);
        // Rows of the domain are the basis vectors b_i, so a point is
        // k = sum_i p_i b_i = B^T p. Fractional coordinates are p = (B^T)^-1 k.
        vector<Vec> domain_T(dimension);
        for (int i = 0; i < dimension; i++) {
            domain_T[i].dimension = dimension;
            for (int j = 0; j < dimension; j++)
                domain_T[i](j) = f.domain[j][i];
        }
        inv_domain = invertMatrix2(domain_T, dimension);
    } else if (domain_absent(f.domain)) {
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
    } else {
        throw std::runtime_error("Domain matrix is smaller than the field dimension ("
                                 + std::to_string(dimension) + ")");
    }

    int rank = f.rank();
    if (rank == 0) {
        // Regular scalar/vector field
        // Check if flat arrays are populated (from HDF5 loading), otherwise use data variant
        if (!f.real_values.empty()) {
            data = transform_data_from_flat(f.real_values, f.imag_values, f.dimension, f.is_complex);
        } else {
            data = transform_data(f.data, f.dimension);
        }
    } else {
        // Indexed field - load into indexed_data structures
        // Check if flat arrays are populated (from HDF5 loading), otherwise use data variant
        if (!f.real_values.empty()) {
            load_indexed_data(f);
        } else {
            load_indexed_data_from_variant(f);
        }
    }
}

void DataEvaluator::load_indexed_data(BaseData& f) {
    // Data layout in flat arrays (w-k ordering):
    // Values stored as: [w0_k0_idx0, w0_k0_idx1, ..., w0_k1_idx0, ..., w1_k0_idx0, ...]
    // Where idx iterates over all tensor indices in row-major order

    // Determine spatial dimensions
    int n_spatial = spatial_count(f);
    int n_w = f.with_w ? f.w_points.size() : 1;
    int rank = f.rank();

    // Calculate total index size (product of all tensor dimensions)
    int idx_size = f.total_index_size();

    if (rank == 1) {
        // Reshape: indexed_data_1d[w][spatial][index] (w-k ordering)
        int d0 = f.inds[0];
        indexed_data_1d.resize(n_w);
        for (int w = 0; w < n_w; w++) {
            indexed_data_1d[w].resize(n_spatial);
            for (int s = 0; s < n_spatial; s++) {
                indexed_data_1d[w][s].resize(d0);
                for (int i = 0; i < d0; i++) {
                    int flat_idx = (w * n_spatial + s) * idx_size + i;
                    float real_val = (flat_idx < f.real_values.size()) ? f.real_values[flat_idx] : 0.0f;
                    float imag_val = (f.is_complex && flat_idx < f.imag_values.size()) ? f.imag_values[flat_idx] : 0.0f;
                    indexed_data_1d[w][s][i] = cfloat(real_val, imag_val);
                }
            }
        }
    } else if (rank == 2) {
        // Reshape: indexed_data_2d[w][spatial][i][j] (w-k ordering)
        int d0 = f.inds[0];
        int d1 = f.inds[1];
        indexed_data_2d.resize(n_w);
        for (int w = 0; w < n_w; w++) {
            indexed_data_2d[w].resize(n_spatial);
            for (int s = 0; s < n_spatial; s++) {
                indexed_data_2d[w][s].resize(d0, vector<cfloat>(d1));
                for (int i = 0; i < d0; i++) {
                    for (int j = 0; j < d1; j++) {
                        int local_idx = i * d1 + j;
                        int flat_idx = (w * n_spatial + s) * idx_size + local_idx;
                        float real_val = (flat_idx < f.real_values.size()) ? f.real_values[flat_idx] : 0.0f;
                        float imag_val = (f.is_complex && flat_idx < f.imag_values.size()) ? f.imag_values[flat_idx] : 0.0f;
                        indexed_data_2d[w][s][i][j] = cfloat(real_val, imag_val);
                    }
                }
            }
        }
    } else if (rank == 3) {
        // Reshape: indexed_data_3d[w][spatial][i][j][k] (w-k ordering)
        int d0 = f.inds[0];
        int d1 = f.inds[1];
        int d2 = f.inds[2];
        indexed_data_3d.resize(n_w);
        for (int w = 0; w < n_w; w++) {
            indexed_data_3d[w].resize(n_spatial);
            for (int s = 0; s < n_spatial; s++) {
                indexed_data_3d[w][s].resize(d0, vector<vector<cfloat>>(d1, vector<cfloat>(d2)));
                for (int i = 0; i < d0; i++) {
                    for (int j = 0; j < d1; j++) {
                        for (int k = 0; k < d2; k++) {
                            int local_idx = (i * d1 + j) * d2 + k;
                            int flat_idx = (w * n_spatial + s) * idx_size + local_idx;
                            float real_val = (flat_idx < f.real_values.size()) ? f.real_values[flat_idx] : 0.0f;
                            float imag_val = (f.is_complex && flat_idx < f.imag_values.size()) ? f.imag_values[flat_idx] : 0.0f;
                            indexed_data_3d[w][s][i][j][k] = cfloat(real_val, imag_val);
                        }
                    }
                }
            }
        }
    } else if (rank == 4) {
        // Reshape: indexed_data_4d[w][spatial][i][j][k][l] (w-k ordering)
        int d0 = f.inds[0];
        int d1 = f.inds[1];
        int d2 = f.inds[2];
        int d3 = f.inds[3];
        indexed_data_4d.resize(n_w);
        for (int w = 0; w < n_w; w++) {
            indexed_data_4d[w].resize(n_spatial);
            for (int s = 0; s < n_spatial; s++) {
                indexed_data_4d[w][s].resize(d0, vector<vector<vector<cfloat>>>(d1,
                    vector<vector<cfloat>>(d2, vector<cfloat>(d3))));
                for (int i = 0; i < d0; i++) {
                    for (int j = 0; j < d1; j++) {
                        for (int k = 0; k < d2; k++) {
                            for (int l = 0; l < d3; l++) {
                                int local_idx = ((i * d1 + j) * d2 + k) * d3 + l;
                                int flat_idx = (w * n_spatial + s) * idx_size + local_idx;
                                float real_val = (flat_idx < f.real_values.size()) ? f.real_values[flat_idx] : 0.0f;
                                float imag_val = (f.is_complex && flat_idx < f.imag_values.size()) ? f.imag_values[flat_idx] : 0.0f;
                                indexed_data_4d[w][s][i][j][k][l] = cfloat(real_val, imag_val);
                            }
                        }
                    }
                }
            }
        }
    }
}

void DataEvaluator::load_indexed_data_from_variant(BaseData& f) {
    // Data layout expected in BaseData (w-k ordering):
    // For rank=1: vector<vector<cfloat>> where outer is w/spatial, inner is index dim
    // For rank=2: vector<vector<vector<cfloat>>> where [w/spatial][i][j]
    // For rank=3: vector<vector<vector<vector<cfloat>>>> where [w/spatial][i][j][k]
    // For rank=4: vector<vector<vector<vector<vector<cfloat>>>>> where [w/spatial][i][j][k][l]

    // Determine spatial dimensions
    int n_spatial = spatial_count(f);
    int n_w = f.with_w ? f.w_points.size() : 1;
    int rank = f.rank();

    if (rank == 1) {
        // Extract 1D indexed data
        if (auto* vec_data = std::get_if<vector<vector<cfloat>>>(&f.data)) {
            // Reshape: indexed_data_1d[w][spatial][index] (w-k ordering)
            indexed_data_1d.resize(n_w);
            for (int w = 0; w < n_w; w++) {
                indexed_data_1d[w].resize(n_spatial);
                for (int s = 0; s < n_spatial; s++) {
                    int flat_idx = w * n_spatial + s;
                    if (flat_idx < vec_data->size()) {
                        indexed_data_1d[w][s] = (*vec_data)[flat_idx];
                    }
                }
            }
        }
    } else if (rank == 2) {
        // Extract 2D indexed data (matrix)
        if (auto* mat_data = std::get_if<vector<vector<vector<cfloat>>>>(&f.data)) {
            // Reshape: indexed_data_2d[w][spatial][i][j] (w-k ordering)
            indexed_data_2d.resize(n_w);
            for (int w = 0; w < n_w; w++) {
                indexed_data_2d[w].resize(n_spatial);
                for (int s = 0; s < n_spatial; s++) {
                    int flat_idx = w * n_spatial + s;
                    if (flat_idx < mat_data->size()) {
                        indexed_data_2d[w][s] = (*mat_data)[flat_idx];
                    }
                }
            }
        }
    } else if (rank == 3) {
        // Extract 3D indexed data
        if (auto* tensor_data = std::get_if<vector<vector<vector<vector<cfloat>>>>>(&f.data)) {
            // Reshape: indexed_data_3d[w][spatial][i][j][k] (w-k ordering)
            indexed_data_3d.resize(n_w);
            for (int w = 0; w < n_w; w++) {
                indexed_data_3d[w].resize(n_spatial);
                for (int s = 0; s < n_spatial; s++) {
                    int flat_idx = w * n_spatial + s;
                    if (flat_idx < tensor_data->size()) {
                        indexed_data_3d[w][s] = (*tensor_data)[flat_idx];
                    }
                }
            }
        }
    } else if (rank == 4) {
        // Extract 4D indexed data
        if (auto* tensor_data = std::get_if<vector<vector<vector<vector<vector<cfloat>>>>>>(&f.data)) {
            // Reshape: indexed_data_4d[w][spatial][i][j][k][l] (w-k ordering)
            indexed_data_4d.resize(n_w);
            for (int w = 0; w < n_w; w++) {
                indexed_data_4d[w].resize(n_spatial);
                for (int s = 0; s < n_spatial; s++) {
                    int flat_idx = w * n_spatial + s;
                    if (flat_idx < tensor_data->size()) {
                        indexed_data_4d[w][s] = (*tensor_data)[flat_idx];
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
        }
        augmented[i][n + i] = 1.0f; // Identity matrix
    }
    // Scale-aware singularity threshold
    float scale = 0.0f;
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            scale = std::max(scale, std::abs(augmented[i][j]));
    float tol = 1e-6f * (scale > 0.0f ? scale : 1.0f);
    // Gaussian Elimination
    for (size_t i = 0; i < n; ++i) {
        // Partial Pivoting
        size_t maxRow = i;
        for (size_t k = i + 1; k < n; ++k) {
            if (std::abs(augmented[k][i]) > std::abs(augmented[maxRow][i])) {
                maxRow = k;
            }
        }
        std::swap(augmented[i], augmented[maxRow]);

        // Check for singular matrix
        if (std::abs(augmented[i][i]) < tol) {
            std::string error_msg = "Matrix is singular and cannot be inverted.\n";
            error_msg += "Matrix after partial pivoting at step " + std::to_string(i) + ":\n";
            for (size_t row = 0; row < n; ++row) {
                error_msg += "  [";
                for (size_t col = 0; col < n; ++col) {
                    error_msg += std::to_string(augmented[row][col]);
                    if (col < n - 1) error_msg += ", ";
                }
                error_msg += "]\n";
            }
            error_msg += "Pivot element [" + std::to_string(i) + "][" + std::to_string(i) + "] = "
                       + std::to_string(augmented[i][i]) + " (too small, threshold = " + std::to_string(tol) + ")";
            throw std::runtime_error(error_msg);
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
        inv[i].dimension = n;
        for (size_t j = 0; j < n; ++j) {
            inv[i](j) = augmented[i][n + j];
        }
    }

    return inv;
}

// ============================================================================
// Grid interpolation
//
// Spatial axes are periodic: an axis with n samples stores sample i at
// fractional coordinate i/n, and interpolation between sample n-1 and the
// next period's sample 0 wraps around. Flat spatial index is row-major over
// mesh (last axis fastest). The frequency axis is non-periodic, varies
// slowest, and must contain the requested w (within a small tolerance).
// ============================================================================

namespace {

struct Stencil {
    int i0 = 0, i1 = 0;
    float r = 0.0f;  // weight of i1
};

Stencil periodic_stencil(float x, int n) {
    Stencil s;
    if (n <= 1) return s;
    float t = x * n;
    float fl = std::floor(t);
    s.r = t - fl;
    long base = (long)fl;
    s.i0 = (int)(((base % n) + n) % n);
    if (s.r >= 1.0f) {  // guard against float rounding
        s.r = 0.0f;
        s.i0 = (s.i0 + 1) % n;
    }
    s.i1 = (s.i0 + 1) % n;
    return s;
}

Stencil frequency_stencil(float w, vector<float>& w_points) {
    Stencil s;
    int n = w_points.size();
    if (n <= 1) return s;
    float w_min = w_points[0], w_max = w_points[n - 1];
    w = sanitize_within_bounds(w, w_min, w_max);
    if (w < w_min || w > w_max)
        throw out_of_range("w=" + std::to_string(w) + " outside frequency range ["
                           + std::to_string(w_min) + ", " + std::to_string(w_max) + "]");
    int i = binary_search(w, w_points);
    if (i < 0) i = 0;
    if (i >= n - 1) i = n - 2;
    s.i0 = i;
    s.i1 = i + 1;
    s.r = (w - w_points[i]) / (w_points[i + 1] - w_points[i]);
    return s;
}

// zero_like / axpy work on scalars, complex<Vec>, and nested vectors of them
inline cfloat zero_like(const cfloat&) { return cfloat(0.0f, 0.0f); }
inline complex<Vec> zero_like(const complex<Vec>& v) {
    Vec z = v.real();
    z.x = z.y = z.z = z.w = 0.0f;
    return complex<Vec>(z, z);
}
template <typename U> vector<U> zero_like(const vector<U>& v) {
    vector<U> out;
    out.reserve(v.size());
    for (auto& e : v) out.push_back(zero_like(e));
    return out;
}

inline void axpy(cfloat& acc, float wt, const cfloat& v) { acc += wt * v; }
inline void axpy(complex<Vec>& acc, float wt, const complex<Vec>& v) { acc = acc + wt * v; }
template <typename U> void axpy(vector<U>& acc, float wt, const vector<U>& v) {
    if (acc.size() != v.size())
        throw std::runtime_error("Inconsistent tensor shapes in field data");
    for (size_t i = 0; i < acc.size(); i++) axpy(acc[i], wt, v[i]);
}

// Multilinear interpolation over the spatial mesh at fractional point p.
// get(spatial_index) returns the stored value at that grid point.
template <typename T, typename Get>
T interpolate_mesh(const vector<int>& mesh, const Vec& p, Get get) {
    int naxes = std::min<int>(mesh.size(), 3);
    if ((int)mesh.size() > 3)
        throw std::runtime_error("Meshes with more than 3 spatial axes are not supported");
    Stencil st[3];
    int stride[3] = {1, 1, 1};
    for (int d = naxes - 1; d >= 0; d--) {
        st[d] = periodic_stencil(p(d), mesh[d]);
        if (d < naxes - 1) stride[d] = stride[d + 1] * mesh[d + 1];
    }

    const T& first = get(0);
    T acc = zero_like(first);
    for (int corner = 0; corner < (1 << naxes); corner++) {
        float weight = 1.0f;
        int idx = 0;
        for (int d = 0; d < naxes; d++) {
            bool upper = (corner >> d) & 1;
            weight *= upper ? st[d].r : 1.0f - st[d].r;
            idx += (upper ? st[d].i1 : st[d].i0) * stride[d];
        }
        if (weight == 0.0f) continue;
        axpy(acc, weight, get(idx));
    }
    return acc;
}

} // namespace

int DataEvaluator::n_spatial() const {
    int n = 1;
    for (int m : mesh) n *= m;
    if (n <= 0) throw std::runtime_error("Field mesh has a zero-length axis");
    return n;
}

Vec DataEvaluator::to_fractional(const Vec& point) const {
    Vec p;
    p.dimension = dimension;
    bool identity = (int)inv_domain.size() < dimension;
    for (int i = 0; i < dimension && i < 3; i++) {
        p(i) = 0;
        for (int j = 0; j < dimension && j < 3; j++) {
            float a = identity ? (i == j ? 1.0f : 0.0f) : inv_domain[i](j);
            p(i) += a * point(j);
        }
        p(i) -= std::floor(p(i));  // fold into first zone [0, 1)
    }
    return p;
}

// src(w_idx, spatial_idx) -> const T&
template <typename T, typename Src>
T DataEvaluator::evaluate(const Vec& p, float w, Src src) const {
    auto at_w = [&](int wi) {
        return interpolate_mesh<T>(mesh, p, [&](int s) -> const T& { return src(wi, s); });
    };
    if (!with_w) return at_w(0);
    Stencil fs = frequency_stencil(w, const_cast<vector<float>&>(w_points));
    T v0 = at_w(fs.i0);
    if (fs.r == 0.0f) return v0;
    T acc = zero_like(v0);
    axpy(acc, 1.0f - fs.r, v0);
    axpy(acc, fs.r, at_w(fs.i1));
    return acc;
}

// Frequency-only evaluation of a k-dependent field: average over all k
template <typename T, typename Src>
T DataEvaluator::evaluate_k_average(float w, Src src) const {
    int nk = n_spatial();
    auto avg_at_w = [&](int wi) {
        T acc = zero_like(src(wi, 0));
        for (int s = 0; s < nk; s++) axpy(acc, 1.0f / nk, src(wi, s));
        return acc;
    };
    if (!with_w) return avg_at_w(0);
    Stencil fs = frequency_stencil(w, const_cast<vector<float>&>(w_points));
    T v0 = avg_at_w(fs.i0);
    if (fs.r == 0.0f) return v0;
    T acc = zero_like(v0);
    axpy(acc, 1.0f - fs.r, v0);
    axpy(acc, fs.r, avg_at_w(fs.i1));
    return acc;
}

static void require(bool ok, const string& msg) {
    if (!ok) throw std::runtime_error(msg);
}

ResultVariant DataEvaluator::operator()(float w) {
    require(!data.empty(), "Cannot evaluate field: no scalar data available");
    int nk = n_spatial();
    int nw = with_w ? w_points.size() : 1;
    require((int)data.size() >= nk * nw, "Field data size (" + to_string(data.size())
            + ") smaller than mesh x frequency points (" + to_string(nk * nw) + ")");
    auto src = [&](int wi, int s) -> const complex<Vec>& { return data[wi * nk + s]; };
    complex<Vec> ans = evaluate_k_average<complex<Vec>>(w, src);
    return convert(ans);
}

ResultVariant DataEvaluator::operator()(Vec point, float w) {
    require(!scattered, "Evaluating fields stored on scattered k-points (as_mesh=false) is not supported");
    require(!data.empty(), "Cannot evaluate field: no scalar data available");
    int nk = n_spatial();
    int nw = with_w ? w_points.size() : 1;
    require((int)data.size() >= nk * nw, "Field data size (" + to_string(data.size())
            + ") smaller than mesh x frequency points (" + to_string(nk * nw) + ")");

    Vec p = to_fractional(point);
    auto src = [&](int wi, int s) -> const complex<Vec>& { return data[wi * nk + s]; };
    complex<Vec> ans = evaluate<complex<Vec>>(p, w, src);
    return convert(ans);
}

static float real_part(const cfloat& v) { return v.real(); }
static vector<float> real_part(const vector<cfloat>& v) {
    vector<float> out(v.size());
    for (size_t i = 0; i < v.size(); i++) out[i] = v[i].real();
    return out;
}
template <typename U> static auto real_part(const vector<vector<U>>& v) {
    vector<decltype(real_part(v[0]))> out;
    out.reserve(v.size());
    for (auto& e : v) out.push_back(real_part(e));
    return out;
}

template <typename T> static ResultVariant as_result(T&& v, bool is_complex) {
    if (is_complex) return ResultVariant(std::forward<T>(v));
    return ResultVariant(real_part(v));
}

ResultVariant DataEvaluator::get_array_impl(const Vec* point, float w) {
    int rank = inds.size();
    if (rank == 0) {
        return point ? (*this)(*point, w) : (*this)(w);
    }
    require(rank <= 4, "Tensor rank > 4 not yet supported");
    require(!scattered, "Evaluating fields stored on scattered k-points (as_mesh=false) is not supported");

    int nk = n_spatial();
    int nw = with_w ? w_points.size() : 1;
    Vec p = point ? to_fractional(*point) : Vec();

    auto run = [&](auto& store) -> ResultVariant {
        using T = std::decay_t<decltype(store[0][0])>;
        require((int)store.size() >= nw, "Indexed field has fewer frequency slices than w_points");
        for (int wi = 0; wi < nw; wi++)
            require((int)store[wi].size() >= nk, "Indexed field has fewer spatial points than mesh");
        auto src = [&](int wi, int s) -> const T& { return store[wi][s]; };
        T v = point ? evaluate<T>(p, w, src) : evaluate_k_average<T>(w, src);
        return as_result(std::move(v), is_complex);
    };

    if ((rank == 1 && indexed_data_1d.empty()) || (rank == 2 && indexed_data_2d.empty()) ||
        (rank == 3 && indexed_data_3d.empty()) || (rank == 4 && indexed_data_4d.empty())) {
        // No data - return an empty array of the right type
        if (rank == 1) return as_result(vector<cfloat>(), is_complex);
        if (rank == 2) return as_result(vector<vector<cfloat>>(), is_complex);
        if (rank == 3) return as_result(vector<vector<vector<cfloat>>>(), is_complex);
        return as_result(vector<vector<vector<vector<cfloat>>>>(), is_complex);
    }
    if (rank == 1) return run(indexed_data_1d);
    if (rank == 2) return run(indexed_data_2d);
    if (rank == 3) return run(indexed_data_3d);
    return run(indexed_data_4d);
}

ResultVariant DataEvaluator::get_array(Vec point, float w) {
    return get_array_impl(&point, w);
}

ResultVariant DataEvaluator::operator()(Vec point, vector<int> indices, float w) {
    int rank = inds.size();

    if (rank == 0) {
        throw std::runtime_error("Indexed operator called on non-indexed field");
    }

    if (indices.size() != rank) {
        throw std::runtime_error("Number of indices doesn't match tensor rank");
    }

    // Validate indices against the actual dimensions (supporting non-uniform tensors)
    for (int i = 0; i < rank; i++) {
        if (indices[i] < 0 || indices[i] >= inds[i]) {
            throw std::runtime_error("Index out of bounds for dimension " + std::to_string(i) +
                                     " (index=" + std::to_string(indices[i]) +
                                     ", max=" + std::to_string(inds[i] - 1) + ")");
        }
    }

    // Evaluate the full tensor and pick out the requested element
    ResultVariant full = get_array(point, w);
    auto pick = [&](auto& arr) -> ResultVariant {
        using A = std::decay_t<decltype(arr)>;
        if constexpr (std::is_same_v<A, vector<cfloat>> || std::is_same_v<A, vector<float>>) {
            return vector<typename A::value_type>{arr[indices[0]]};
        } else if constexpr (std::is_same_v<A, vector<vector<cfloat>>> || std::is_same_v<A, vector<vector<float>>>) {
            return vector<typename A::value_type::value_type>{arr[indices[0]][indices[1]]};
        } else if constexpr (std::is_same_v<A, vector<vector<vector<cfloat>>>> || std::is_same_v<A, vector<vector<vector<float>>>>) {
            return vector<typename A::value_type::value_type::value_type>{arr[indices[0]][indices[1]][indices[2]]};
        } else if constexpr (std::is_same_v<A, vector<vector<vector<vector<cfloat>>>>> || std::is_same_v<A, vector<vector<vector<vector<float>>>>>) {
            return vector<typename A::value_type::value_type::value_type::value_type>{arr[indices[0]][indices[1]][indices[2]][indices[3]]};
        } else {
            throw std::runtime_error("Indexed evaluation returned a non-array result");
        }
    };
    return std::visit(pick, full);
}

ResultVariant DataEvaluator::get_array(float w) {
    return get_array_impl(nullptr, w);
}
