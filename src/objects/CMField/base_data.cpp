// base_field.cpp
#include "base_data.hpp"
#include <H5Cpp.h>
#include <stdexcept>
#include <iostream>
#include <cmath>

using namespace H5;

void read_vector(vector<int> &vec, DataSet &ds) {
    try {
        DataSpace space = ds.getSpace();
        hsize_t size;
        space.getSimpleExtentDims(&size);
        vec.resize(size);
        ds.read(vec.data(), PredType::NATIVE_INT);
        space.close();
        ds.close();
    } catch (...) {
        // Legacy format: use empty inds for scalar
        vec = {};
    }
}

void read_vector(vector<float> &vec, DataSet &ds) {
    try {
        DataSpace space = ds.getSpace();
        hsize_t size;
        space.getSimpleExtentDims(&size);
        vec.resize(size);
        ds.read(vec.data(), PredType::NATIVE_FLOAT);
        space.close();
        ds.close();
    } catch (...) {
        // Legacy format: use empty inds for scalar
        vec = {};
    }
}
void load_metadata(BaseData &field, H5File &file) {
    // Read bools as ints to avoid size mismatch issues
    int temp_is_complex, temp_is_vector, temp_with_k, temp_with_w, temp_as_mesh;
    file.openDataSet("/is_complex").read(&temp_is_complex, PredType::NATIVE_INT);
    file.openDataSet("/is_vector").read(&temp_is_vector, PredType::NATIVE_INT);
    file.openDataSet("/with_k").read(&temp_with_k, PredType::NATIVE_INT);
    file.openDataSet("/with_w").read(&temp_with_w, PredType::NATIVE_INT);

    field.is_complex = temp_is_complex;
    field.is_vector = temp_is_vector;
    field.with_k = temp_with_k;
    field.with_w = temp_with_w;

    file.openDataSet("/dimension").read(&field.dimension, PredType::NATIVE_INT);

    file.openDataSet("/as_mesh").read(&temp_as_mesh, PredType::NATIVE_INT);
    field.as_mesh = temp_as_mesh;

    // Read centered (with fallback for legacy files)
    try {
        int temp_centered;
        file.openDataSet("/centered").read(&temp_centered, PredType::NATIVE_INT);
        field.centered = temp_centered;
    } catch (...) {
        field.centered = true;  // Default for legacy files
    }

    // Read inds (tensor indices)
    DataSet inds_ds = file.openDataSet("/inds");
    read_vector(field.inds, inds_ds);
}

// ============================================================================
// Helper functions for writing (defined before use)
// ============================================================================

void write_metadata(H5File& file, bool is_complex, bool is_vector, bool with_k, bool with_w, bool as_mesh, bool centered, int dimension, const vector<int>& inds, const vector<int>& mesh) {
    auto write_scalar = [&](const std::string& name, int value) {
        DataSpace scalar_space(H5S_SCALAR);
        DataSet ds = file.createDataSet(name, PredType::NATIVE_INT, scalar_space);
        ds.write(&value, PredType::NATIVE_INT);
    };

    write_scalar("/is_complex", is_complex);
    write_scalar("/is_vector", is_vector);
    write_scalar("/with_k", with_k);
    write_scalar("/with_w", with_w);
    write_scalar("/dimension", dimension);
    write_scalar("/as_mesh", as_mesh);
    write_scalar("/centered", centered);

    // Write inds array
    if (!inds.empty()) {
        hsize_t dims[1] = { inds.size() };
        DataSpace space(1, dims);
        DataSet ds = file.createDataSet("/inds", PredType::NATIVE_INT, space);
        ds.write(inds.data(), PredType::NATIVE_INT);
    } else {
        hsize_t dims[1] = { 0 };
        DataSpace space(1, dims);
        file.createDataSet("/inds", PredType::NATIVE_INT, space);
    }

    // Write mesh array
    if (!mesh.empty()) {
        hsize_t dims[1] = { mesh.size() };
        DataSpace space(1, dims);
        DataSet ds = file.createDataSet("/mesh", PredType::NATIVE_INT, space);
        ds.write(mesh.data(), PredType::NATIVE_INT);
    }
}

void write_points(H5File& file, const std::vector<std::vector<float>>& points) {
    hsize_t dims[2] = { points.size(), points[0].size() };
    DataSpace space(2, dims);
    DataSet ds = file.createDataSet("/points", PredType::NATIVE_FLOAT, space);

    // flatten 2D into 1D buffer
    std::vector<float> flat;
    flat.reserve(points.size() * points[0].size());
    for (auto const& row : points) {
        flat.insert(flat.end(), row.begin(), row.end());
    }
    ds.write(flat.data(), PredType::NATIVE_FLOAT);
}

void write_domain(H5File& file, const std::vector<std::vector<float>>& domain) {
    hsize_t dims[2] = { domain.size(), domain[0].size() };
    DataSpace space(2, dims);
    DataSet ds = file.createDataSet("/domain", PredType::NATIVE_FLOAT, space);

    // flatten 2D into 1D buffer
    std::vector<float> flat;
    flat.reserve(domain.size() * domain[0].size());
    for (auto const& row : domain) {
        flat.insert(flat.end(), row.begin(), row.end());
    }
    ds.write(flat.data(), PredType::NATIVE_FLOAT);
}

void write_w_points(H5File& file, const std::vector<float>& w_points) {
    hsize_t dims[1] = { w_points.size() };
    DataSpace space(1, dims);
    DataSet ds = file.createDataSet("/w_points", PredType::NATIVE_FLOAT, space);
    ds.write(w_points.data(), PredType::NATIVE_FLOAT);
}

// Flatten functions for complex types
void flatten(const vector<vector<vector<vector<cfloat>>>>& tensor, vector<float>& real_flat, vector<float>& imag_flat, bool is_complex) {
    for (auto& vol : tensor) {
        for (auto& plane : vol) {
            for (auto& row : plane) {
                for (auto& v : row) {
                    real_flat.push_back(v.real());
                    if (is_complex) imag_flat.push_back(v.imag());
                }
            }
        }
    }
}

void flatten(const vector<vector<vector<cfloat>>>& tensor, vector<float>& real_flat, vector<float>& imag_flat, bool is_complex) {
    for (auto& plane : tensor) {
        for (auto& row : plane) {
            for (auto& v : row) {
                real_flat.push_back(v.real());
                if (is_complex) imag_flat.push_back(v.imag());
            }
        }
    }
}

void flatten(const vector<vector<cfloat>>& tensor, vector<float>& real_flat, vector<float>& imag_flat, bool is_complex) {
    for (auto& row : tensor) {
        for (auto& v : row) {
            real_flat.push_back(v.real());
            if (is_complex) imag_flat.push_back(v.imag());
        }
    }
}

void flatten(const vector<cfloat>& tensor, vector<float>& real_flat, vector<float>& imag_flat, bool is_complex) {
    for (auto& v : tensor) {
        real_flat.push_back(v.real());
        if (is_complex) imag_flat.push_back(v.imag());
    }
}

void flatten(const vector<vector<vector<vector<vector<cfloat>>>>>& tensor, vector<float>& real_flat, vector<float>& imag_flat, bool is_complex) {
    for (auto& hyper : tensor) {
        for (auto& vol : hyper) {
            for (auto& plane : vol) {
                for (auto& row : plane) {
                    for (auto& v : row) {
                        real_flat.push_back(v.real());
                        if (is_complex) imag_flat.push_back(v.imag());
                    }
                }
            }
        }
    }
}

// Flatten functions for float types (no .real() needed)
void flatten(const vector<vector<vector<vector<float>>>>& tensor, vector<float>& real_flat) {
    for (auto& vol : tensor) {
        for (auto& plane : vol) {
            for (auto& row : plane) {
                for (auto& v : row) {
                    real_flat.push_back(v);
                }
            }
        }
    }
}

void flatten(const vector<vector<vector<float>>>& tensor, vector<float>& real_flat) {
    for (auto& plane : tensor) {
        for (auto& row : plane) {
            for (auto& v : row) {
                real_flat.push_back(v);
            }
        }
    }
}

void flatten(const vector<vector<float>>& tensor, vector<float>& real_flat) {
    for (auto& row : tensor) {
        for (auto& v : row) {
            real_flat.push_back(v);
        }
    }
}

// ============================================================================
// Reading functions
// ============================================================================

void store_domain(DataSet& ds_domain, BaseData& field) {
    DataSpace space_domain = ds_domain.getSpace();
    int rank = space_domain.getSimpleExtentNdims();
    if (rank != 2) {
        printf("Error: Domain rank = %d. Exiting\n", rank);
        exit(1);
    }
    hsize_t dims[2];
    space_domain.getSimpleExtentDims(dims);
    // Read into flat buffer first, then populate 2D structure
    std::vector<float> domain_flat(dims[0] * dims[1]);
    ds_domain.read(domain_flat.data(), PredType::NATIVE_FLOAT);
    field.domain.assign(dims[0], std::vector<float>(dims[1]));
    for (size_t i = 0; i < dims[0]; i++) {
        for (size_t j = 0; j < dims[1]; j++) {
            field.domain[i][j] = domain_flat[i * dims[1] + j];
            if (std::isnan(field.domain[i][j]) || std::isinf(field.domain[i][j])) {
                throw std::runtime_error("Domain matrix contains NaN or Inf values");
            }
        }
    }
    space_domain.close();
    ds_domain.close();
}

void load_k_points(DataSet& ds_points, BaseData& field) {
    DataSpace space = ds_points.getSpace();
    int rank = space.getSimpleExtentNdims();
    if (rank != 2) {
        printf("Error: Points rank = %d. Exiting\n", rank);
        exit(1);
    }
    hsize_t dims[2];
    space.getSimpleExtentDims(dims);
    // Read into flat buffer first, then populate 2D structure
    std::vector<float> points_flat(dims[0] * dims[1]);
    ds_points.read(points_flat.data(), PredType::NATIVE_FLOAT);
    field.points.assign(dims[0], std::vector<float>(dims[1]));
    for (size_t i = 0; i < dims[0]; i++) {
        for (size_t j = 0; j < dims[1]; j++) {
            field.points[i][j] = points_flat[i * dims[1] + j];
        }
    }
    space.close();
    ds_points.close();
}

BaseData load_data_from_hdf5(const std::string& filename) {
    BaseData field;
    //Check if file exists
    if (FILE *file = fopen(filename.c_str(), "r")) {
        fclose(file);
    } else {
        throw std::runtime_error("File not found: " + filename);
    }
    H5File file(filename, H5F_ACC_RDONLY);

    // -- Metadata --
    load_metadata(field, file);

    // -- Domain (optional) --
    if (field.with_k && field.as_mesh) {
        try {
            DataSet ds_domain = file.openDataSet("/domain");
            store_domain(ds_domain, field);
        } catch (...) {
            printf("Warning: /domain dataset is missing or invalid. Treating as empty domain.\n");
        }
    }

    // -- Mesh (optional) --
    if (field.as_mesh) {
        try {
            DataSet ds = file.openDataSet("/mesh");
            read_vector(field.mesh, ds);
        } catch (...) {
            printf("Warning: /mesh dataset is missing or invalid. Treating as empty mesh.\n");
        }
    }

    // -- Points (k-point data when as_mesh = false) --
    if (!field.as_mesh && field.with_k) {
        try {
            DataSet ds = file.openDataSet("/points");
            load_k_points(ds, field);
        } catch (...) {
            printf("Warning: /points dataset is missing or invalid. Treating as empty points.\n");
        }
    }

    // -- w_points (optional) --
    if (field.with_w) {
        try {
            DataSet ds = file.openDataSet("/w_points");
            read_vector(field.w_points, ds);
        } catch (...) {
            printf("Warning: /w_points dataset is missing or invalid. Treating as empty w_points.\n");
        }
    }

    // --- Real part ---
    H5::DataSet real_ds = file.openDataSet("/values/real");
    read_vector(field.real_values, real_ds);

    // --- Imag part (if complex) ---
    if (field.is_complex) {
        H5::DataSet imag_ds = file.openDataSet("/values/imag");
        read_vector(field.imag_values, imag_ds);  // Store in field, not local var
    }

    // Explicitly close the file to release locks immediately
    file.close();
    return field;
}

// Load with ordering specification (ordering ignored - w-k ordering is assumed)
BaseData load_data_from_hdf5(const std::string& filename, const std::string& ordering) {
    // Ordering parameter is deprecated - w-k ordering is always used
    return load_data_from_hdf5(filename);
}

// ============================================================================
// Saving functions
// ============================================================================

void save_data_to_hdf5(const std::string& filename,
                        bool is_complex,
                        bool is_vector,
                        bool with_k,
                        bool with_w,
                        bool as_mesh,
                        bool centered,
                        const std::vector<int>& inds,
                        const std::vector<int>& mesh,
                        const std::vector<std::vector<float>>& domain,
                        int dimension,
                        const std::vector<float>& w_points,
                        const std::vector<std::vector<float>>& points,
                        const vector<vector<vector<vector<cfloat>>>>& data) {

    H5File file(filename, H5F_ACC_TRUNC);
    write_metadata(file, is_complex, is_vector, with_k, with_w, as_mesh, centered, dimension, inds, mesh);

    if (!points.empty()) write_points(file, points);
    if (!domain.empty()) write_domain(file, domain);
    if (!w_points.empty()) write_w_points(file, w_points);

    std::vector<float> real_flat, imag_flat;
    flatten(data, real_flat, imag_flat, is_complex);

    // -- Write datasets under /values --
    hsize_t dims[1] = { real_flat.size() };
    DataSpace space(1, dims);

    Group values_group = file.createGroup("/values");
    // Real part
    DataSet ds_real = values_group.createDataSet("real", PredType::NATIVE_FLOAT, space);
    ds_real.write(real_flat.data(), PredType::NATIVE_FLOAT);

    // Imag part
    if (is_complex) {
        DataSet ds_imag = values_group.createDataSet("imag", PredType::NATIVE_FLOAT, space);
        ds_imag.write(imag_flat.data(), PredType::NATIVE_FLOAT);
    }
}

void save_data_to_hdf5(const std::string& filename,
                        bool is_complex,
                        bool is_vector,
                        bool with_k,
                        bool with_w,
                        bool as_mesh,
                        bool centered,
                        const std::vector<int>& inds,
                        const std::vector<int>& mesh,
                        const std::vector<std::vector<float>>& domain,
                        int dimension,
                        const std::vector<float>& w_points,
                        const std::vector<std::vector<float>>& points,
                        const vector<vector<vector<cfloat>>>& data) {

    H5File file(filename, H5F_ACC_TRUNC);
    write_metadata(file, is_complex, is_vector, with_k, with_w, as_mesh, centered, dimension, inds, mesh);

    if (!points.empty()) write_points(file, points);
    if (!domain.empty()) write_domain(file, domain);
    if (!w_points.empty()) write_w_points(file, w_points);

    std::vector<float> real_flat, imag_flat;
    flatten(data, real_flat, imag_flat, is_complex);

    // -- Write datasets under /values --
    hsize_t dims[1] = { real_flat.size() };
    DataSpace space(1, dims);

    Group values_group = file.createGroup("/values");
    // Real part
    DataSet ds_real = values_group.createDataSet("real", PredType::NATIVE_FLOAT, space);
    ds_real.write(real_flat.data(), PredType::NATIVE_FLOAT);

    // Imag part
    if (is_complex) {
        DataSet ds_imag = values_group.createDataSet("imag", PredType::NATIVE_FLOAT, space);
        ds_imag.write(imag_flat.data(), PredType::NATIVE_FLOAT);
    }
}

void save_data_to_hdf5(const std::string& filename,
                        bool is_complex,
                        bool is_vector,
                        bool with_k,
                        bool with_w,
                        bool as_mesh,
                        bool centered,
                        const std::vector<int>& inds,
                        const std::vector<int>& mesh,
                        const std::vector<std::vector<float>>& domain,
                        int dimension,
                        const std::vector<float>& w_points,
                        const std::vector<std::vector<float>>& points,
                        const vector<vector<cfloat>>& data) {

    H5File file(filename, H5F_ACC_TRUNC);
    write_metadata(file, is_complex, is_vector, with_k, with_w, as_mesh, centered, dimension, inds, mesh);

    if (!points.empty()) write_points(file, points);
    if (!domain.empty()) write_domain(file, domain);
    if (!w_points.empty()) write_w_points(file, w_points);

    std::vector<float> real_flat, imag_flat;
    flatten(data, real_flat, imag_flat, is_complex);

    // -- Write datasets under /values --
    hsize_t dims[1] = { real_flat.size() };
    DataSpace space(1, dims);

    Group values_group = file.createGroup("/values");
    // Real part
    DataSet ds_real = values_group.createDataSet("real", PredType::NATIVE_FLOAT, space);
    ds_real.write(real_flat.data(), PredType::NATIVE_FLOAT);

    // Imag part
    if (is_complex) {
        DataSet ds_imag = values_group.createDataSet("imag", PredType::NATIVE_FLOAT, space);
        ds_imag.write(imag_flat.data(), PredType::NATIVE_FLOAT);
    }
}

void save_data_to_hdf5(const std::string& filename,
                        bool is_complex,
                        bool is_vector,
                        bool with_k,
                        bool with_w,
                        bool as_mesh,
                        bool centered,
                        const std::vector<int>& inds,
                        const std::vector<int>& mesh,
                        const std::vector<std::vector<float>>& domain,
                        int dimension,
                        const std::vector<float>& w_points,
                        const std::vector<std::vector<float>>& points,
                        const vector<cfloat>& data) {

    H5File file(filename, H5F_ACC_TRUNC);
    write_metadata(file, is_complex, is_vector, with_k, with_w, as_mesh, centered, dimension, inds, mesh);

    if (!points.empty()) write_points(file, points);
    if (!domain.empty()) write_domain(file, domain);
    if (!w_points.empty()) write_w_points(file, w_points);

    std::vector<float> real_flat, imag_flat;
    flatten(data, real_flat, imag_flat, is_complex);

    // -- Write datasets under /values --
    hsize_t dims[1] = { real_flat.size() };
    DataSpace space(1, dims);

    Group values_group = file.createGroup("/values");
    // Real part
    DataSet ds_real = values_group.createDataSet("real", PredType::NATIVE_FLOAT, space);
    ds_real.write(real_flat.data(), PredType::NATIVE_FLOAT);

    // Imag part
    if (is_complex) {
        DataSet ds_imag = values_group.createDataSet("imag", PredType::NATIVE_FLOAT, space);
        ds_imag.write(imag_flat.data(), PredType::NATIVE_FLOAT);
    }
}

void save_data_to_hdf5(const std::string& filename,
                        bool is_complex,
                        bool is_vector,
                        bool with_k,
                        bool with_w,
                        bool as_mesh,
                        bool centered,
                        const std::vector<int>& inds,
                        const std::vector<int>& mesh,
                        const std::vector<std::vector<float>>& domain,
                        int dimension,
                        const std::vector<float>& w_points,
                        const std::vector<std::vector<float>>& points,
                        const vector<vector<vector<vector<float>>>>& data) {

    H5File file(filename, H5F_ACC_TRUNC);
    write_metadata(file, is_complex, is_vector, with_k, with_w, as_mesh, centered, dimension, inds, mesh);

    if (!points.empty()) write_points(file, points);
    if (!domain.empty()) write_domain(file, domain);
    if (!w_points.empty()) write_w_points(file, w_points);

    std::vector<float> real_flat;
    flatten(data, real_flat);

    // -- Write datasets under /values --
    hsize_t dims[1] = { real_flat.size() };
    DataSpace space(1, dims);

    Group values_group = file.createGroup("/values");
    // Real part
    DataSet ds_real = values_group.createDataSet("real", PredType::NATIVE_FLOAT, space);
    ds_real.write(real_flat.data(), PredType::NATIVE_FLOAT);
}

void save_data_to_hdf5(const std::string& filename,
                        bool is_complex,
                        bool is_vector,
                        bool with_k,
                        bool with_w,
                        bool as_mesh,
                        bool centered,
                        const std::vector<int>& inds,
                        const std::vector<int>& mesh,
                        const std::vector<std::vector<float>>& domain,
                        int dimension,
                        const std::vector<float>& w_points,
                        const std::vector<std::vector<float>>& points,
                        const vector<vector<vector<float>>>& data) {

    H5File file(filename, H5F_ACC_TRUNC);
    write_metadata(file, is_complex, is_vector, with_k, with_w, as_mesh, centered, dimension, inds, mesh);

    if (!points.empty()) write_points(file, points);
    if (!domain.empty()) write_domain(file, domain);
    if (!w_points.empty()) write_w_points(file, w_points);

    std::vector<float> real_flat;
    flatten(data, real_flat);

    // -- Write datasets under /values --
    hsize_t dims[1] = { real_flat.size() };
    DataSpace space(1, dims);

    Group values_group = file.createGroup("/values");
    // Real part
    DataSet ds_real = values_group.createDataSet("real", PredType::NATIVE_FLOAT, space);
    ds_real.write(real_flat.data(), PredType::NATIVE_FLOAT);
}

void save_data_to_hdf5(const std::string& filename,
                        bool is_complex,
                        bool is_vector,
                        bool with_k,
                        bool with_w,
                        bool as_mesh,
                        bool centered,
                        const std::vector<int>& inds,
                        const std::vector<int>& mesh,
                        const std::vector<std::vector<float>>& domain,
                        int dimension,
                        const std::vector<float>& w_points,
                        const std::vector<std::vector<float>>& points,
                        const vector<vector<float>>& data) {

    H5File file(filename, H5F_ACC_TRUNC);
    write_metadata(file, is_complex, is_vector, with_k, with_w, as_mesh, centered, dimension, inds, mesh);

    if (!points.empty()) write_points(file, points);
    if (!domain.empty()) write_domain(file, domain);
    if (!w_points.empty()) write_w_points(file, w_points);

    std::vector<float> real_flat;
    flatten(data, real_flat);

    // -- Write datasets under /values --
    hsize_t dims[1] = { real_flat.size() };
    DataSpace space(1, dims);

    Group values_group = file.createGroup("/values");
    // Real part
    DataSet ds_real = values_group.createDataSet("real", PredType::NATIVE_FLOAT, space);
    ds_real.write(real_flat.data(), PredType::NATIVE_FLOAT);
}

void save_data_to_hdf5(const std::string& filename,
                        bool is_complex,
                        bool is_vector,
                        bool with_k,
                        bool with_w,
                        bool as_mesh,
                        bool centered,
                        const std::vector<int>& inds,
                        const std::vector<int>& mesh,
                        const std::vector<std::vector<float>>& domain,
                        int dimension,
                        const std::vector<float>& w_points,
                        const std::vector<std::vector<float>>& points,
                        const vector<float>& data) {

    H5File file(filename, H5F_ACC_TRUNC);
    write_metadata(file, is_complex, is_vector, with_k, with_w, as_mesh, centered, dimension, inds, mesh);

    if (!points.empty()) write_points(file, points);
    if (!domain.empty()) write_domain(file, domain);
    if (!w_points.empty()) write_w_points(file, w_points);

    std::vector<float> real_flat = data;

    // -- Write datasets under /values --
    hsize_t dims[1] = { real_flat.size() };
    DataSpace space(1, dims);

    Group values_group = file.createGroup("/values");
    // Real part
    DataSet ds_real = values_group.createDataSet("real", PredType::NATIVE_FLOAT, space);
    ds_real.write(real_flat.data(), PredType::NATIVE_FLOAT);
}

void save_data_to_hdf5(const std::string& filename,
                        bool is_complex,
                        bool is_vector,
                        bool with_k,
                        bool with_w,
                        bool as_mesh,
                        bool centered,
                        const std::vector<int>& inds,
                        const std::vector<int>& mesh,
                        const std::vector<std::vector<float>>& domain,
                        int dimension,
                        const std::vector<float>& w_points,
                        const std::vector<std::vector<float>>& points,
                        const vector<vector<vector<vector<vector<cfloat>>>>>& data) {

    H5File file(filename, H5F_ACC_TRUNC);
    write_metadata(file, is_complex, is_vector, with_k, with_w, as_mesh, centered, dimension, inds, mesh);

    if (!points.empty()) write_points(file, points);
    if (!domain.empty()) write_domain(file, domain);
    if (!w_points.empty()) write_w_points(file, w_points);

    std::vector<float> real_flat, imag_flat;
    flatten(data, real_flat, imag_flat, is_complex);

    // -- Write datasets under /values --
    hsize_t dims[1] = { real_flat.size() };
    DataSpace space(1, dims);

    Group values_group = file.createGroup("/values");
    // Real part
    DataSet ds_real = values_group.createDataSet("real", PredType::NATIVE_FLOAT, space);
    ds_real.write(real_flat.data(), PredType::NATIVE_FLOAT);

    // Imag part
    if (is_complex) {
        DataSet ds_imag = values_group.createDataSet("imag", PredType::NATIVE_FLOAT, space);
        ds_imag.write(imag_flat.data(), PredType::NATIVE_FLOAT);
    }
}

void save_data_to_hdf5(BaseData& field, const std::string& filename) {
    std::visit([&](auto& data) {
        save_data_to_hdf5(filename, field.is_complex, field.is_vector,
                          field.with_k, field.with_w, field.as_mesh,
                          field.centered, field.inds, field.mesh, field.domain,
                          field.dimension, field.w_points, field.points, data);
    }, field.data);
}

// ============================================================================
// Convenience save_data wrappers
// ============================================================================

void save_data(string filename, vector<vector<vector<vector<cfloat>>>>& data, vector<int> inds, vector<int> mesh, vector<vector<float>> domain, vector<float> w_points, const vector<vector<float>>& points, bool centered) {
    bool with_k = points.size() > 0 || mesh.size() > 0;
    save_data_to_hdf5(filename, true, false, with_k, w_points.size() > 0, mesh.size() > 0, centered, inds, mesh, domain, mesh.size(), w_points, points, data);
}

void save_data(string filename, vector<vector<vector<cfloat>>>& data, vector<int> inds, vector<int> mesh, vector<vector<float>> domain, vector<float> w_points, const vector<vector<float>>& points, bool centered) {
    bool with_k = points.size() > 0 || mesh.size() > 0;
    save_data_to_hdf5(filename, true, false, with_k, w_points.size() > 0, mesh.size() > 0, centered, inds, mesh, domain, mesh.size(), w_points, points, data);
}

void save_data(string filename, vector<vector<cfloat>>& data, vector<int> inds, vector<int> mesh, vector<vector<float>> domain, vector<float> w_points, const vector<vector<float>>& points, bool centered) {
    bool with_k = points.size() > 0 || mesh.size() > 0;
    save_data_to_hdf5(filename, true, false, true, w_points.size() > 0, mesh.size() > 0, centered, inds, mesh, domain, mesh.size(), w_points, points, data);
}

void save_data(string filename, vector<cfloat>& data, vector<int> inds, vector<int> mesh, vector<vector<float>> domain, vector<float> w_points, const vector<vector<float>>& points, bool centered) {
    bool with_k = points.size() > 0 || mesh.size() > 0;
    save_data_to_hdf5(filename, true, false, with_k, w_points.size() > 0, mesh.size() > 0, centered, inds, mesh, domain, mesh.size(), w_points, points, data);
}

void save_data(string filename, vector<vector<vector<vector<float>>>>& data, vector<int> inds, vector<int> mesh, vector<vector<float>> domain, vector<float> w_points, const vector<vector<float>>& points, bool centered) {
    bool with_k = points.size() > 0 || mesh.size() > 0;
    save_data_to_hdf5(filename, false, false, with_k, w_points.size() > 0, mesh.size() > 0, centered, inds, mesh, domain, mesh.size(), w_points, points, data);
}

void save_data(string filename, vector<vector<vector<float>>>& data, vector<int> inds, vector<int> mesh, vector<vector<float>> domain, vector<float> w_points, const vector<vector<float>>& points, bool centered) {
    bool with_k = points.size() > 0 || mesh.size() > 0;
    save_data_to_hdf5(filename, false, false, with_k, w_points.size() > 0, mesh.size() > 0, centered, inds, mesh, domain, mesh.size(), w_points, points, data);
}

void save_data(string filename, vector<vector<float>>& data, vector<int> inds, vector<int> mesh, vector<vector<float>> domain, vector<float> w_points, const vector<vector<float>>& points, bool centered) {
    bool with_k = points.size() > 0 || mesh.size() > 0;
    save_data_to_hdf5(filename, false, false, true, w_points.size() > 0, mesh.size() > 0, centered, inds, mesh, domain, mesh.size(), w_points, points, data);
}

void save_data(string filename, vector<float>& data, vector<int> inds, vector<int> mesh, vector<vector<float>> domain, vector<float> w_points, const vector<vector<float>>& points, bool centered) {
    bool with_k = points.size() > 0 || mesh.size() > 0;
    save_data_to_hdf5(filename, false, false, with_k, w_points.size() > 0, mesh.size() > 0, centered, inds, mesh, domain, mesh.size(), w_points, points, data);
}

void save_data(string filename, vector<float>& data, vector<float> w_points, bool centered) {
    vector<int> tmp;
    vector<vector<float>> tmp2;
    save_data_to_hdf5(filename, false, false, false, w_points.size() > 0, false, centered, tmp, tmp, tmp2, 0, w_points, tmp2, data);
}
