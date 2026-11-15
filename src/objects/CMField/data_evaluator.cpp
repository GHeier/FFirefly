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

    n_points = spatial_count(f);
    block = inds.empty() ? f.vec_len() : f.total_index_size();

    vector<float> re, im;
    f.flat_values(re, im);
    values.resize(re.size());
    for (size_t i = 0; i < re.size(); i++)
        values[i] = cfloat(re[i], i < im.size() ? im[i] : 0.0f);
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

} // namespace

int DataEvaluator::n_spatial() const {
    if (n_points <= 0) throw std::runtime_error("Field mesh has a zero-length axis");
    return n_points;
}

int DataEvaluator::n_freq() const {
    return with_w ? std::max<int>(w_points.size(), 1) : 1;
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

// Grid points and weights for one evaluation: multilinear over the spatial
// mesh at `point` (or a uniform average over all k when point is null),
// then linear in frequency.
DataEvaluator::Taps DataEvaluator::taps(const Vec* point, float w) const {
    int nk = n_spatial();
    Taps spatial;
    if (point) {
        int naxes = std::min<int>(mesh.size(), 3);
        if ((int)mesh.size() > 3)
            throw std::runtime_error("Meshes with more than 3 spatial axes are not supported");
        Vec p = to_fractional(*point);
        Stencil st[3];
        long stride[3] = {1, 1, 1};
        for (int d = naxes - 1; d >= 0; d--) {
            st[d] = periodic_stencil(p(d), mesh[d]);
            if (d < naxes - 1) stride[d] = stride[d + 1] * mesh[d + 1];
        }
        for (int corner = 0; corner < (1 << naxes); corner++) {
            float weight = 1.0f;
            long idx = 0;
            for (int d = 0; d < naxes; d++) {
                bool upper = (corner >> d) & 1;
                weight *= upper ? st[d].r : 1.0f - st[d].r;
                idx += (upper ? st[d].i1 : st[d].i0) * stride[d];
            }
            if (weight != 0.0f) spatial.push_back({idx, weight});
        }
    } else {
        spatial.reserve(nk);
        for (int s = 0; s < nk; s++) spatial.push_back({s, 1.0f / nk});
    }
    if (!with_w) return spatial;

    Stencil fs = frequency_stencil(w, const_cast<vector<float>&>(w_points));
    Taps out;
    out.reserve(2 * spatial.size());
    for (auto& [s, wt] : spatial) {
        out.push_back({(long)fs.i0 * nk + s, wt * (1.0f - fs.r)});
        if (fs.r != 0.0f) out.push_back({(long)fs.i1 * nk + s, wt * fs.r});
    }
    return out;
}

// Weighted sum of `count` consecutive elements starting at `offset` within each grid point's block
vector<cfloat> DataEvaluator::combine(const Taps& t, int offset, int count) const {
    long needed = (long)n_freq() * n_spatial() * block;
    if ((long)values.size() < needed)
        throw std::runtime_error("Field data size (" + to_string(values.size())
                                 + ") smaller than frequency x mesh x values per point ("
                                 + to_string(needed) + ")");
    vector<cfloat> out(count, cfloat(0.0f, 0.0f));
    for (auto& [g, wt] : t) {
        const cfloat* src = values.data() + g * block + offset;
        for (int e = 0; e < count; e++) out[e] += wt * src[e];
    }
    return out;
}

static void require(bool ok, const string& msg) {
    if (!ok) throw std::runtime_error(msg);
}

// Scalar fields return float or cfloat. Complex vector fields return
// complex<Vec>; real vector fields keep the historical behaviour of returning
// the first component as a float.
ResultVariant DataEvaluator::scalar_result(const vector<cfloat>& v) const {
    if (is_complex && is_vector) {
        Vec re, im;
        re.dimension = im.dimension = dimension;
        for (int i = 0; i < (int)v.size() && i < 3; i++) {
            re(i) = v[i].real();
            im(i) = v[i].imag();
        }
        return complex<Vec>(re, im);
    }
    if (is_complex) return v[0];
    return v[0].real();
}

namespace {

template <typename T> T element(const cfloat& c);
template <> cfloat element<cfloat>(const cfloat& c) { return c; }
template <> float element<float>(const cfloat& c) { return c.real(); }

template <typename T, int R> struct Nested { using type = vector<typename Nested<T, R - 1>::type>; };
template <typename T> struct Nested<T, 0> { using type = T; };

// Row-major reshape of a flat block into nested vectors shaped by inds
template <typename T, int R>
typename Nested<T, R>::type reshape(const cfloat*& p, const vector<int>& inds, int axis) {
    if constexpr (R == 0) {
        return element<T>(*p++);
    } else {
        typename Nested<T, R>::type out;
        out.reserve(inds[axis]);
        for (int i = 0; i < inds[axis]; i++) out.push_back(reshape<T, R - 1>(p, inds, axis + 1));
        return out;
    }
}

template <typename T, int R>
ResultVariant shaped(const vector<cfloat>& v, const vector<int>& inds) {
    if (v.empty()) return typename Nested<T, R>::type{};
    const cfloat* p = v.data();
    return reshape<T, R>(p, inds, 0);
}

template <typename T>
ResultVariant shaped_rank(const vector<cfloat>& v, const vector<int>& inds) {
    switch (inds.size()) {
        case 1: return shaped<T, 1>(v, inds);
        case 2: return shaped<T, 2>(v, inds);
        case 3: return shaped<T, 3>(v, inds);
        case 4: return shaped<T, 4>(v, inds);
    }
    throw std::runtime_error("Tensor rank > 4 not yet supported");
}

} // namespace

ResultVariant DataEvaluator::tensor_result(const vector<cfloat>& v) const {
    return is_complex ? shaped_rank<cfloat>(v, inds) : shaped_rank<float>(v, inds);
}

ResultVariant DataEvaluator::operator()(float w) {
    require(inds.empty(), "Use get_array() to evaluate a tensor field");
    require(!values.empty(), "Cannot evaluate field: no scalar data available");
    return scalar_result(combine(taps(nullptr, w), 0, block));
}

ResultVariant DataEvaluator::operator()(Vec point, float w) {
    require(inds.empty(), "Use get_array() to evaluate a tensor field");
    require(!scattered, "Evaluating fields stored on scattered k-points (as_mesh=false) is not supported");
    require(!values.empty(), "Cannot evaluate field: no scalar data available");
    return scalar_result(combine(taps(&point, w), 0, block));
}

ResultVariant DataEvaluator::get_array(Vec point, float w) {
    if (inds.empty()) return (*this)(point, w);
    require(inds.size() <= 4, "Tensor rank > 4 not yet supported");
    require(!scattered, "Evaluating fields stored on scattered k-points (as_mesh=false) is not supported");
    if (values.empty()) return tensor_result({});
    return tensor_result(combine(taps(&point, w), 0, block));
}

ResultVariant DataEvaluator::get_array(float w) {
    if (inds.empty()) return (*this)(w);
    require(inds.size() <= 4, "Tensor rank > 4 not yet supported");
    if (values.empty()) return tensor_result({});
    return tensor_result(combine(taps(nullptr, w), 0, block));
}

// A single tensor element, returned as a one-element vector
ResultVariant DataEvaluator::operator()(Vec point, vector<int> indices, float w) {
    int rank = inds.size();
    require(rank > 0, "Indexed operator called on non-indexed field");
    require((int)indices.size() == rank, "Number of indices doesn't match tensor rank");
    require(!scattered, "Evaluating fields stored on scattered k-points (as_mesh=false) is not supported");

    int offset = 0;
    for (int i = 0; i < rank; i++) {
        if (indices[i] < 0 || indices[i] >= inds[i]) {
            throw std::runtime_error("Index out of bounds for dimension " + std::to_string(i) +
                                     " (index=" + std::to_string(indices[i]) +
                                     ", max=" + std::to_string(inds[i] - 1) + ")");
        }
        offset = offset * inds[i] + indices[i];
    }
    require(!values.empty(), "Cannot evaluate field: no data available");
    cfloat v = combine(taps(&point, w), offset, 1)[0];
    if (is_complex) return vector<cfloat>{v};
    return vector<float>{v.real()};
}
