#pragma once
#include <vector>
#include <complex>
#include <string>
#include <variant>
#include <map>
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

    // Flat value arrays, w-major: [(w_idx * nk + k_idx) * block + element].
    // Filled when loading from HDF5, and by flatten_variant() for fields
    // built from `data`. These are what evaluation and saving read.
    vector<float> real_values;
    vector<float> imag_values;

    // Structured data variant (programmatic creation; empty after loading)
    DataVariant data;

    // Free-form metadata, stored as string attributes on the HDF5 root group.
    // Conventional keys: quantity, units, frequency_axis (real | matsubara),
    // temperature, created_by. See docs in CLAUDE.md (HDF5 Data Format).
    map<string, string> attrs;

    // HDF5 layout version. Files written before versioning are read as 1.
    static constexpr int FORMAT_VERSION = 2;
    int format_version = FORMAT_VERSION;

    // Values as flat real/imag arrays (imag is empty for real fields). The
    // variant `data` is used when it holds values, otherwise real_values/imag_values.
    void flat_values(vector<float>& re, vector<float>& im) const;

    // Copy `data` into real_values/imag_values if they are still empty
    void flatten_variant();

    // Value of attrs[key], or `fallback` when it is not set
    string attr(const string& key, const string& fallback = "") const {
        auto it = attrs.find(key);
        return it == attrs.end() ? fallback : it->second;
    }

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
            // Version 1 files with with_w but no w_points stored nw in mesh[0]
            int start_idx = (format_version < 2 && with_w && w_points.empty()) ? 1 : 0;
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
        // Version 1 files stored nw in mesh[0] (matching nk(), which skips it)
        if (format_version < 2 && as_mesh && !mesh.empty()) return mesh[0];
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
// Spatial dimension implied by a mesh (trailing singleton axes dropped, e.g.
// {24, 24, 1} -> 2) or, without a mesh, by the length of the k-points.
int infer_dimension(const vector<int>& mesh, const vector<vector<float>>& points);

// Save overloads
void save_data_to_hdf5(BaseData& data, const std::string& filename);

// Read or add metadata attributes on an existing field file (see BaseData::attrs)
map<string, string> read_field_attributes(const std::string& filename);
void write_field_attributes(const std::string& filename, const map<string, string>& attrs);
void save_data(string filename, vector<vector<vector<vector<cfloat>>>>& data, vector<int> inds = {}, vector<int> mesh = {}, vector<vector<float>> domain = {{}}, vector<float> w_points = {}, const vector<vector<float>>& points = ep(), bool centered = true);
void save_data(string filename, vector<vector<vector<cfloat>>>& data, vector<int> inds = {}, vector<int> mesh = {}, vector<vector<float>> domain = {{}}, vector<float> w_points = {}, const vector<vector<float>>& points = ep(), bool centered = true);
void save_data(string filename, vector<vector<cfloat>>& data, vector<int> inds = {}, vector<int> mesh = {}, vector<vector<float>> domain = {{}}, vector<float> w_points = {}, const vector<vector<float>>& points = ep(), bool centered = true);
void save_data(string filename, vector<cfloat>& data, vector<int> inds = {}, vector<int> mesh = {}, vector<vector<float>> domain = {{}}, vector<float> w_points = {}, const vector<vector<float>>& points = ep(), bool centered = true);

void save_data(string filename, vector<vector<vector<vector<float>>>>& data, vector<int> inds = {}, vector<int> mesh = {}, vector<vector<float>> domain = {{}}, vector<float> w_points = {}, const vector<vector<float>>& points = ep(), bool centered = true);
void save_data(string filename, vector<vector<vector<float>>>& data, vector<int> inds = {}, vector<int> mesh = {}, vector<vector<float>> domain = {{}}, vector<float> w_points = {}, const vector<vector<float>>& points = ep(), bool centered = true);
void save_data(string filename, vector<vector<float>>& data, vector<int> inds = {}, vector<int> mesh = {}, vector<vector<float>> domain = {{}}, vector<float> w_points = {}, const vector<vector<float>>& points = ep(), bool centered = true);
void save_data(string filename, vector<float>& data, vector<int> inds = {}, vector<int> mesh = {}, vector<vector<float>> domain = {{}}, vector<float> w_points = {}, const vector<vector<float>>& points = ep(), bool centered = true);
void save_data(string filename, vector<float>& data, vector<float> w_points, bool centered = true);

