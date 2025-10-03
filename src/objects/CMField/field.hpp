#pragma once

#include "base_data.hpp"
#include "data_evaluator.hpp"
#include "../vec.hpp"
#include <vector>
#include <complex>
#include <cmath>

using namespace std;
using cfloat = complex<float>;

class Field {
public:
    BaseData data;

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
            // Get the domain extent for this dimension
            float extent = 0.0;
            for (int j = 0; j < dimension; j++) {
                extent += data.domain[i][j] * data.domain[i][j];
            }
            extent = sqrt(extent);

            // Wrap to [-extent/2, extent/2]
            while (wrapped(i) > extent / 2.0) wrapped(i) -= extent;
            while (wrapped(i) < -extent / 2.0) wrapped(i) += extent;
        }
        return wrapped;
    }

    // Helper to initialize shift vectors and evaluator
    void initialize() {
        // Infer dimension from domain
        dimension = data.domain.empty() ? 1 : data.domain.size();
        data.dimension = dimension;

        // Create shift vectors for centering
        shift_vectors.resize(1);
        shift_vectors[0].dimension = dimension;

        for (int i = 0; i < dimension; i++) {
            // Calculate shift to center the domain
            float shift = 0.0;
            for (int j = 0; j < dimension; j++) {
                shift += data.domain[i][j] * 0.5;
            }
            shift_vectors[0](i) = shift;
        }

        // Initialize evaluator
        evaluator = DataEvaluator(data);
    }

public:
    // Constructor with default values
    Field(const BaseData::DataVariant& data_variant,
          bool is_complex = false,
          bool is_vector = false,
          const vector<int>& mesh = {},
          const vector<vector<float>>& domain = {},
          const vector<float>& w_points = {})
    {
        data.data = data_variant;
        data.is_complex = is_complex;
        data.is_vector = is_vector;
        data.mesh = mesh;
        data.domain = domain;
        data.w_points = w_points;

        // Infer with_w based on w_points
        data.with_w = !w_points.empty();

        // Infer with_k based on mesh
        data.with_k = !mesh.empty();

        // Infer as_mesh
        data.as_mesh = !mesh.empty();

        initialize();
    }

    // Constructor from file
    Field(const string& filename) {
        data = load_data_from_hdf5(filename);
        initialize();
    }

    // Operator for frequency-only evaluation
    ResultVariant operator()(float w) {
        return evaluator(w);
    }

    // Operator for spatial evaluation with optional frequency
    ResultVariant operator()(Vec point, float w = 0) {
        // Apply periodic boundary conditions
        Vec periodic_point = apply_periodic_bc(point);

        // Shift to DataEvaluator's coordinate system (which starts at 0)
        Vec shifted_point = periodic_point;
        for (int i = 0; i < dimension; i++) {
            shifted_point(i) = periodic_point(i) + shift_vectors[0](i);
        }

        return evaluator(shifted_point, w);
    }

    // Copy assignment operator
    Field& operator=(const Field& other) {
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
