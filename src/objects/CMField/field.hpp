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
    vector<Vec> shift_vectors;  // Store shift for each dimension
    int dimension;

    // Helper function to shift point to centered coordinates
    Vec shift_to_centered(Vec point) const {
        Vec shifted = point;
        for (int i = 0; i < dimension; i++) {
            shifted(i) = point(i) - shift_vectors[0](i);
        }
        return shifted;
    }

    // Helper function to apply periodic boundary conditions
    Vec apply_periodic_bc(Vec point) const {
        Vec wrapped = point;
        for (int i = 0; i < dimension; i++) {
            // Get the domain extent for this dimension (i-th diagonal element)
            float extent = data.domain[i][i];

            // Wrap to [-extent/2, extent/2]
            while (wrapped(i) > extent / 2.0) wrapped(i) -= extent;
            while (wrapped(i) < -extent / 2.0) wrapped(i) += extent;
        }
        return wrapped;
    }

    // Helper to initialize shift vectors and evaluator
    void initialize() {
        // Prioritize stored dimension, then infer from domain or mesh
        if (data.dimension > 0) {
            // Use stored dimension (e.g., from HDF5 file)
            dimension = data.dimension;
        } else if (!data.domain.empty()) {
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

        // Create shift vectors for centering
        shift_vectors.resize(1);
        shift_vectors[0].dimension = dimension;

        // Only calculate shift if centered mode is enabled
        if (centered && !data.domain.empty()) {
            for (int i = 0; i < dimension; i++) {
                // Calculate shift to center the domain
                // Use only the diagonal element (i-th component of i-th basis vector)
                shift_vectors[0](i) = data.domain[i][i] * 0.5;
            }
        } else {
            // No centering - zero shift
            for (int i = 0; i < dimension; i++) {
                shift_vectors[0](i) = 0.0;
            }
        }

        // Initialize evaluator
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
    FieldImpl(const string& filename, bool centered_coords = true) {
        data = load_data_from_hdf5(filename);
        centered = centered_coords;
        initialize();
    }

    // Operator for frequency-only evaluation
    ResultVariant operator()(float w) {
        return evaluator(w);
    }

    // Operator for spatial evaluation with optional frequency
    ResultVariant operator()(Vec point, float w = 0) {
        if (!data.with_k) return evaluator(w);
        // Apply periodic boundary conditions
        Vec periodic_point = apply_periodic_bc(point);

        // Shift to DataEvaluator's coordinate system (which starts at 0)
        Vec shifted_point = periodic_point;
        for (int i = 0; i < dimension; i++) {
            shifted_point(i) = periodic_point(i) + shift_vectors[0](i);
        }

        return evaluator(shifted_point, w);
    }

    // Operator for indexed spatial evaluation (e.g., matrix elements H_ab(k))
    ResultVariant operator()(Vec point, vector<int> indices, float w = 0) {
        if (!data.with_k) return evaluator(w);
        // Apply periodic boundary conditions
        Vec periodic_point = apply_periodic_bc(point);

        // Shift to DataEvaluator's coordinate system (which starts at 0)
        Vec shifted_point = periodic_point;
        for (int i = 0; i < dimension; i++) {
            shifted_point(i) = periodic_point(i) + shift_vectors[0](i);
        }

        return evaluator(shifted_point, indices, w);
    }

    // Get full array at point (e.g., H(k) returns full matrix)
    ResultVariant get_array(Vec point, float w = 0) {
        // Apply periodic boundary conditions
        Vec periodic_point = apply_periodic_bc(point);

        // Shift to DataEvaluator's coordinate system (which starts at 0)
        Vec shifted_point = periodic_point;
        for (int i = 0; i < dimension; i++) {
            shifted_point(i) = periodic_point(i) + shift_vectors[0](i);
        }

        return evaluator.get_array(shifted_point, w);
    }

    // Copy assignment operator
    FieldImpl& operator=(const FieldImpl& other) {
        if (this != &other) {
            data = other.data;
            shift_vectors = other.shift_vectors;
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
