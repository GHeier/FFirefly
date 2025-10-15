#pragma once
#include <vector>
#include <complex>
#include <string>
#include <variant>
#include <cmath>
#include <stdexcept>
#include "../vec.hpp"

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

    int n_indices = 0;
    int dim_indices = 1;
    vector<int> mesh;
    int dimension = 1;
    vector<vector<float>> domain;

    vector<float> w_points;
    vector<vector<float>> points;  // k-point data storage when as_mesh = false

    using DataVariant = variant<
        vector<cfloat>,
        vector<vector<cfloat>>,
        vector<vector<vector<cfloat>>>,  // For n_indices arrays at each point
        vector<vector<vector<vector<cfloat>>>>  // For nested multi-dimensional indices
    >;

    DataVariant data;

    int total_index_size() const { return pow(dim_indices, n_indices); }
    int nk() const {
        if (!with_k) return 1;
        if (as_mesh) {
            int total = 1;
            for (int m : mesh) total *= m;
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

// Save overloads
void save_data_to_hdf5(BaseData& data, const std::string& filename);
void save_data(string filename, BaseData::DataVariant& data, bool is_complex = false, vector<int> mesh = {}, vector<vector<float>> domain = {{}}, vector<float> w_points = {}, int n_indices = 0, int dim_indices = 0);

void save_data_to_hdf5(const std::string& filename, bool is_complex, bool is_vector, bool is_matrix, bool with_k, bool with_w, bool as_mesh, int n_indices, int dim_indices, vector<int> &mesh, vector<vector<float>> &domain, int dimension, vector<float> &w_points, vector<vector<float>> &points, const BaseData::DataVariant& data);
