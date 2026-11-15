#pragma once

#include "base_data.hpp"
#include "data_evaluator.hpp"
#include "src/objects/vec.hpp"
#include <vector>
#include <complex>
#include <cmath>
#include <iostream>

using namespace std;
using cfloat = complex<float>;

// Internal field implementation class
class FieldImpl {
public:
    BaseData data;
    bool centered;  // Whether to apply coordinate centering

private:
    DataEvaluator evaluator;
    Vec shift;      // Cartesian offset from centered coordinates to the grid origin
    int dimension;

    static bool has_domain(const BaseData& d) {
        for (auto& row : d.domain)
            if (!row.empty()) return true;
        return false;
    }

    // Map a (possibly centered) Cartesian point into the evaluator's frame.
    // Periodicity is handled by the evaluator in fractional coordinates.
    Vec to_grid_frame(Vec point) const {
        Vec shifted = point;
        for (int i = 0; i < dimension && i < 3; i++) {
            shifted(i) = point(i) + shift(i);
        }
        return shifted;
    }

    // Helper to initialize shift vector and evaluator
    void initialize() {
        // Prioritize stored dimension, then infer from domain or mesh
        if (data.dimension > 0) {
            // Use stored dimension (e.g., from HDF5 file)
            dimension = data.dimension;
        } else if (has_domain(data)) {
            // Infer from domain matrix size
            dimension = data.domain.size();
            data.dimension = dimension;
        } else if (!data.mesh.empty()) {
            // Infer from mesh, counting only dimensions > 1
            dimension = 0;
            for (int m : data.mesh) {
                if (m > 1) dimension++;
            }
            if (dimension == 0) dimension = 1;  // At least 1D
            data.dimension = dimension;
        } else {
            dimension = 1;
            data.dimension = dimension;
        }
        data.centered = centered;

        // Centered coordinates put the grid origin at -(b_1 + ... + b_d)/2,
        // i.e. fractional coordinate -1/2 along every basis vector.
        shift = Vec();
        shift.dimension = dimension;
        if (centered) {
            Vec half;
            half.dimension = dimension;
            for (int i = 0; i < dimension && i < 3; i++) half(i) = 0.5f;
            if (has_domain(data)) {
                shift = frac_to_cart(data.domain, half);
                shift.dimension = dimension;
            } else {
                shift = half;
            }
        }

        // Keep the flat arrays filled so get_data() works for fields built in code
        data.flatten_variant();
        evaluator = DataEvaluator(data);
    }

public:
    // Constructor with default values
    FieldImpl(const BaseData::DataVariant& data_variant,
          bool is_complex = false,
          bool is_vector = false,
          const vector<int>& mesh = {},
          const vector<vector<float>>& domain = {},
          const vector<float>& w_points = {},
          const vector<int>& inds = {},
          bool centered_coords = true)
    {
        data.data = data_variant;
        data.is_complex = is_complex;
        data.is_vector = is_vector;
        data.mesh = mesh;
        data.domain = domain;
        data.w_points = w_points;
        data.inds = inds;
        data.centered = centered_coords;
        centered = centered_coords;

        // Infer with_w based on w_points
        data.with_w = !w_points.empty();

        // Infer with_k based on mesh
        data.with_k = !mesh.empty();

        // Infer as_mesh
        data.as_mesh = !mesh.empty();

        initialize();
    }

    // Constructor from file
    FieldImpl(const string& filename) {
        data = load_data_from_hdf5(filename);
        centered = data.centered;  // Use the value saved in the file
        initialize();
    }

    // Operator for frequency-only evaluation
    ResultVariant operator()(float w) {
        return evaluator(w);
    }

    // Operator for spatial evaluation with optional frequency
    ResultVariant operator()(Vec point, float w = 0) {
        if (!data.with_k) return evaluator(w);
        return evaluator(to_grid_frame(point), w);
    }

    // Operator for indexed spatial evaluation (e.g., matrix elements H_ab(k))
    ResultVariant operator()(Vec point, vector<int> indices, float w = 0) {
        return evaluator(to_grid_frame(point), indices, w);
    }

    // Get full array at point (e.g., H(k) returns full matrix)
    ResultVariant get_array(Vec point, float w = 0) {
        return evaluator.get_array(to_grid_frame(point), w);
    }

    // Full array as a function of frequency only (k-averaged if k-dependent)
    ResultVariant get_array(float w) {
        return evaluator.get_array(w);
    }

    // Copy assignment operator
    FieldImpl& operator=(const FieldImpl& other) {
        if (this != &other) {
            data = other.data;
            centered = other.centered;
            shift = other.shift;
            dimension = other.dimension;
            evaluator = DataEvaluator(data);
        }
        return *this;
    }

    // Save to file
    void save(const string& filename) {
        save_data_to_hdf5(data, filename);
    }

    // Access to underlying BaseData
    BaseData& get_data() {
        return data;
    }

    const BaseData& get_data() const {
        return data;
    }
};
