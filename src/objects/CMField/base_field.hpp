#pragma once
#include <vector>
#include <complex>
#include <string>
#include <variant>
#include <cmath>
#include <stdexcept>

using cfloat = std::complex<float>;
using namespace std;

struct BaseField {
    bool is_complex = false;
    bool is_vector = false;
    bool with_k = false;
    bool with_w = false;
    bool as_mesh = false;

    int n_indices = 0;       
    int dim_indices = 1;     
    vector<int> mesh;   
    int dimension = 1;       
    vector<float> domain;
    vector<float> inv_domain;

    vector<float> w_points;   

    using DataVariant = variant<
        vector<cfloat>,         
        vector<vector<cfloat>>  
    >;

    DataVariant data;

    int total_index_size() const { return pow(dim_indices, n_indices); }
    int nk() const { return with_k ? mesh[0] : 1; }
    int nw() const { return with_w ? (with_k ? mesh[1] : mesh[0]) : 1; }
    int vec_len() const { return is_vector ? dimension : 1; }
    int total_size() const { return total_index_size() * nk() * nw() * vec_len(); }

    template <typename T>
    T& get() {
        if (auto* ptr = std::get_if<T>(&data)) return *ptr;
        throw std::runtime_error("BaseField: type mismatch");
    }

    template <typename T>
    const T& get() const {
        if (auto* ptr = std::get_if<T>(&data)) return *ptr;
        throw std::runtime_error("BaseField: type mismatch");
    }
};

// Load
BaseField load_field_from_hdf5(const std::string& filename);

// Save overloads
void save_field_to_hdf5(BaseField& field, const std::string& filename);

void save_field_to_hdf5(const std::string& filename, bool is_complex, bool is_vector, bool with_k, bool with_w, bool as_mesh, int n_indices, int dim_indices, vector<int> &mesh, int dimension, vector<float> &w_points, const BaseField::DataVariant& data);
