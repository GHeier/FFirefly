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
    // Data variant type for storing structured data (used for saving)
    using DataVariant = std::variant<
        vector<cfloat>,
        vector<vector<cfloat>>,
        vector<vector<vector<cfloat>>>,
        vector<vector<vector<vector<cfloat>>>>,
        vector<vector<vector<vector<vector<cfloat>>>>>
    >;

    bool is_complex = false;
    bool is_vector = false;
    bool with_k = false;
    bool with_w = false;
    bool as_mesh = false;
    bool centered = true;  // Whether coordinates are centered (default true for backwards compatibility)

    // Tensor indices: inds[i] = size of i-th tensor dimension
    //   Single-band 4-vertex: inds = {1, 1, 1, 1}
    vector<int> inds;

    vector<int> mesh;      // Spatial mesh dimensions (k-points): e.g. {60, 60} for 60x60 grid
    int dimension = 0;     // Spatial dimension (1D, 2D, or 3D); 0 = needs inference
    vector<vector<float>> domain;  // Real-space domain (lattice vectors)

    vector<float> w_points;  // Frequency points (if with_w==true)
    vector<vector<float>> points;  // k-point data storage when as_mesh = false

    // Flat value arrays (populated when loading from HDF5)
    vector<float> real_values;
    vector<float> imag_values;

    // Structured data variant (used for saving, may be empty after loading)
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

//void load_metadata(BaseData &field, H5File &file);

//void read_vector(vector<int> &vec, DataSet &ds);
//void read_vector(vector<float> &vec, DataSet &ds);

//void store_domain(DataSet& ds_domain, BaseData& field);
//void load_k_points(DataSet& ds_points, BaseData& field);
// Load
BaseData load_data_from_hdf5(const std::string& filename);
BaseData load_data_from_hdf5(const std::string& filename, const std::string& ordering);  // ordering: "k-w" or "w-k"
inline const std::vector<std::vector<float>>& ep() {
    static const std::vector<std::vector<float>> e;
    return e;
}
// Save overloads
void save_data_to_hdf5(BaseData& data, const std::string& filename);
void save_data(string filename, vector<vector<vector<vector<cfloat>>>>& data, vector<int> inds = {}, vector<int> mesh = {}, vector<vector<float>> domain = {{}}, vector<float> w_points = {}, const vector<vector<float>>& points = ep(), bool centered = true);
void save_data(string filename, vector<vector<vector<cfloat>>>& data, vector<int> inds = {}, vector<int> mesh = {}, vector<vector<float>> domain = {{}}, vector<float> w_points = {}, const vector<vector<float>>& points = ep(), bool centered = true);
void save_data(string filename, vector<vector<cfloat>>& data, vector<int> inds = {}, vector<int> mesh = {}, vector<vector<float>> domain = {{}}, vector<float> w_points = {}, const vector<vector<float>>& points = ep(), bool centered = true);
void save_data(string filename, vector<cfloat>& data, vector<int> inds = {}, vector<int> mesh = {}, vector<vector<float>> domain = {{}}, vector<float> w_points = {}, const vector<vector<float>>& points = ep(), bool centered = true);

void save_data(string filename, vector<vector<vector<vector<float>>>>& data, vector<int> inds = {}, vector<int> mesh = {}, vector<vector<float>> domain = {{}}, vector<float> w_points = {}, const vector<vector<float>>& points = ep(), bool centered = true);
void save_data(string filename, vector<vector<vector<float>>>& data, vector<int> inds = {}, vector<int> mesh = {}, vector<vector<float>> domain = {{}}, vector<float> w_points = {}, const vector<vector<float>>& points = ep(), bool centered = true);
void save_data(string filename, vector<vector<float>>& data, vector<int> inds = {}, vector<int> mesh = {}, vector<vector<float>> domain = {{}}, vector<float> w_points = {}, const vector<vector<float>>& points = ep(), bool centered = true);
void save_data(string filename, vector<float>& data, vector<int> inds = {}, vector<int> mesh = {}, vector<vector<float>> domain = {{}}, vector<float> w_points = {}, const vector<vector<float>>& points = ep(), bool centered = true);
void save_data(string filename, vector<float>& data, vector<float> w_points, bool centered = true);

//void save_data_to_hdf5(std::string& filename, bool is_complex, bool is_vector, bool with_k, bool with_w, bool as_mesh, const vector<int>& inds, vector<int> &mesh, vector<vector<float>> &domain, int dimension, vector<float> &w_points, vector<vector<float>> &points, const BaseData::DataVariant& data);
