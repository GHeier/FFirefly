#pragma once
#include <vector>
#include <complex>
#include <string>
#include <variant>
#include <cmath>
#include <stdexcept>
#include "src/objects/vec.hpp"

using cfloat = std::complex<float>;
using namespace std;

struct BaseData {
public:
    bool is_complex = false;
    bool is_vector = false;
    bool is_matrix = false;
    bool with_k = false;
    bool with_w = false;
    bool as_mesh = false;

    // Tensor indices: inds[i] = size of i-th tensor dimension
    // Examples:
    //   Scalar: inds = {}
    //   3-vector: inds = {3}
    //   3x3 matrix: inds = {3, 3}
    //   2x3 matrix: inds = {2, 3}
    //   Single-band 4-vertex: inds = {1, 1, 1, 1}
    vector<int> inds;

    vector<int> mesh;      // Spatial mesh dimensions (k-points): e.g. {60, 60} for 60x60 grid
    int dimension = 0;     // Spatial dimension (1D, 2D, or 3D); 0 = needs inference
    vector<vector<float>> domain;  // Real-space domain (lattice vectors)

    vector<float> w_points;  // Frequency points (if with_w==true)
    vector<vector<float>> points;  // k-point data storage when as_mesh = false

    using DataVariant = variant<
        vector<cfloat>,                                // rank=0: scalars (inds={})
        vector<vector<cfloat>>,                        // rank=1: vectors (inds={n})
        vector<vector<vector<cfloat>>>,                // rank=2: matrices (inds={m,n})
        vector<vector<vector<vector<cfloat>>>>,        // rank=3: 3D tensors (inds={l,m,n})
        vector<vector<vector<vector<vector<cfloat>>>>> // rank=4: 4D tensors (inds={i,j,k,l})
    >;

    DataVariant data;

    // Calculate total number of tensor elements
    int total_index_size() const {
        if (inds.empty()) return 1;  // Scalar
        int total = 1;
        for (int d : inds) total *= d;
        return total;
    }

    // Get tensor rank (number of indices)
    int rank() const { return inds.size(); }
    int nk() const {
        if (!with_k) return 1;
        if (as_mesh) {
            int total = 1;
            // Skip first element if with_w AND w_points is empty (legacy: mesh[0] is nw in that case)
            // But if w_points is provided separately, don't skip - mesh contains only spatial dims
            int start_idx = (with_w && w_points.empty()) ? 1 : 0;
            for (int i = start_idx; i < mesh.size(); i++) {
                total *= mesh[i];
            }
            return total;
        }
        return points.size();
    }
    int nw() const {
        if (!with_w) return 1;
        // If w_points is provided, use its size
        if (!w_points.empty()) return w_points.size();
        // Otherwise use mesh dimensions (legacy behavior)
        if (as_mesh) return with_k ? mesh[1] : mesh[0];
        return 1;
    }
    int vec_len() const { return is_vector ? dimension : 1; }
    int total_size() const { return total_index_size() * nk() * nw() * vec_len(); }

    template <typename T>
    T& get() {
        if (auto* ptr = std::get_if<T>(&data)) return *ptr;
        throw std::runtime_error("BaseData: type mismatch");
    }

    template <typename T>
    const T& get() const {
        if (auto* ptr = std::get_if<T>(&data)) return *ptr;
        throw std::runtime_error("BaseData: type mismatch");
    }
};

// Load
BaseData load_data_from_hdf5(const std::string& filename);
BaseData load_data_from_hdf5(const std::string& filename, const std::string& ordering);  // ordering: "k-w" or "w-k"

// Save overloads
void save_data_to_hdf5(BaseData& data, const std::string& filename);
void save_data_to_hdf5(BaseData& data, const std::string& filename, const std::string& ordering);  // ordering: "k-w" or "w-k"
void save_data(string filename, BaseData::DataVariant& data, bool is_complex = false, vector<int> mesh = {}, vector<vector<float>> domain = {{}}, vector<float> w_points = {}, const vector<int>& inds = {});
void save_data_with_points(string filename, BaseData::DataVariant& data, bool is_complex, vector<int> mesh, vector<vector<float>> domain, vector<float> w_points, const vector<int>& inds, vector<vector<float>>& points);

void save_data_to_hdf5(const std::string& filename, bool is_complex, bool is_vector, bool is_matrix, bool with_k, bool with_w, bool as_mesh, const vector<int>& inds, vector<int> &mesh, vector<vector<float>> &domain, int dimension, vector<float> &w_points, vector<vector<float>> &points, const BaseData::DataVariant& data);
