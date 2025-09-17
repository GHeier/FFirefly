// base_field.cpp
#include "base_field.hpp"
#include <H5Cpp.h>
#include <stdexcept>
#include <iostream>

using namespace H5;

BaseField load_field_from_hdf5(const std::string& filename) {
    BaseField field;
    H5File file(filename, H5F_ACC_RDONLY);

    // -- Metadata --
    file.openDataSet("/is_complex").read(&field.is_complex, PredType::NATIVE_INT);
    file.openDataSet("/is_vector").read(&field.is_vector, PredType::NATIVE_INT);
    file.openDataSet("/with_k").read(&field.with_k, PredType::NATIVE_INT);
    file.openDataSet("/with_w").read(&field.with_w, PredType::NATIVE_INT);
    file.openDataSet("/n_indices").read(&field.n_indices, PredType::NATIVE_INT);
    file.openDataSet("/dim_indices").read(&field.dim_indices, PredType::NATIVE_INT);
    file.openDataSet("/dimension").read(&field.dimension, PredType::NATIVE_INT);
    file.openDataSet("/as_mesh").read(&field.as_mesh, PredType::NATIVE_INT);

    DataSet ds_domain = file.openDataSet("/domain");
    DataSpace space = ds.getSpace();
    // -- Mesh --
    if (field.as_mesh) {
        DataSet ds = file.openDataSet("/mesh");
        DataSpace space = ds.getSpace();
        hsize_t dims[1];
        space.getSimpleExtentDims(dims);
        field.mesh.resize(dims[0]);
        ds.read(field.mesh.data(), PredType::NATIVE_INT);
    }

    // -- w_points (optional) --
    if (field.with_w) {
        try {
            DataSet ds = file.openDataSet("/w_points");
            DataSpace space = ds.getSpace();
            hsize_t n_w = 0;
            space.getSimpleExtentDims(&n_w);
            field.w_points.resize(n_w);
            ds.read(field.w_points.data(), PredType::NATIVE_FLOAT);
        } catch (...) {
            // ignore if missing
        }
    }

    // -- Compute sizes --
    int total_indices = pow(field.dim_indices, field.n_indices);
    int nk = field.nk();
    int nw = field.nw();
    int vec_len = field.vec_len();
    int total_elements = total_indices * nk * nw * vec_len;

    // --- Real part ---
    H5::DataSet real_ds = file.openDataSet("/values/real");
    H5::DataSpace real_space = real_ds.getSpace();
    hsize_t dims[1];
    real_space.getSimpleExtentDims(dims);  // number of elements

    std::vector<float> real_flat(dims[0]); // resize to correct size
    real_ds.read(real_flat.data(), H5::PredType::NATIVE_FLOAT);

    // --- Imag part (if complex) ---
    std::vector<float> imag_flat;
    if (field.is_complex) {
        H5::DataSet imag_ds = file.openDataSet("/values/imag");
        H5::DataSpace imag_space = imag_ds.getSpace();
        hsize_t dims_imag[1];
        imag_space.getSimpleExtentDims(dims_imag);

        imag_flat.resize(dims_imag[0]);
        imag_ds.read(imag_flat.data(), H5::PredType::NATIVE_FLOAT);
    }
    // -- Populate variant --
    if (vec_len == 1) {
        // Matrix of scalars (flattened as vector<cfloat>)
        vector<cfloat> flat(total_elements);
        for (int i = 0; i < total_elements; ++i) {
            flat[i] = cfloat(real_flat[i], field.is_complex ? imag_flat[i] : 0.0f);
        }
        field.data = flat;
    }
    else {
        // Matrix of vectors (2D: nk × vec_len)
        vector<vector<cfloat>> mat(nk, vector<cfloat>(vec_len));
        for (int i = 0; i < nk; ++i) {
            for (int v = 0; v < vec_len; ++v) {
                int idx = i * vec_len + v;
                mat[i][v] = cfloat(real_flat[idx], field.is_complex ? imag_flat[idx] : 0.0f);
            }
        }
        field.data = mat;
    }

    return field;
}

void save_field_to_hdf5(BaseField& field, const std::string& filename) {
    save_field_to_hdf5(filename, field.is_complex, field.is_vector, field.with_k, field.with_w, field.as_mesh, field.n_indices, field.dim_indices, field.mesh, field.dimension, field.w_points, field.data);
}

void save_field_to_hdf5(const std::string& filename, bool is_complex, bool is_vector, bool with_k, bool with_w, bool as_mesh, int n_indices, int dim_indices, vector<int> &mesh, int dimension, vector<float> &w_points, const BaseField::DataVariant& data) {
    H5File file(filename, H5F_ACC_TRUNC);

    // -- Write metadata scalars --
    auto write_scalar = [&](const string& name, int value) {
        DataSpace scalar_space(H5S_SCALAR);
        DataSet ds = file.createDataSet(name, PredType::NATIVE_INT, scalar_space);
        ds.write(&value, PredType::NATIVE_INT);
    };

    write_scalar("/is_complex", is_complex);
    write_scalar("/is_vector", is_vector);
    write_scalar("/with_k", with_k);
    write_scalar("/with_w", with_w);
    write_scalar("/n_indices", n_indices);
    write_scalar("/dim_indices", dim_indices);
    write_scalar("/dimension", dimension);
    write_scalar("/as_mesh", as_mesh);

    // -- Mesh --
    if (!mesh.empty()) {
        hsize_t dims[1] = { mesh.size() };
        DataSpace space(1, dims);
        DataSet ds = file.createDataSet("/mesh", PredType::NATIVE_INT, space);
        ds.write(mesh.data(), PredType::NATIVE_INT);
    }

    // -- w_points --
    if (!w_points.empty()) {
        hsize_t dims[1] = { w_points.size() };
        DataSpace space(1, dims);
        DataSet ds = file.createDataSet("/w_points", PredType::NATIVE_FLOAT, space);
        ds.write(w_points.data(), PredType::NATIVE_FLOAT);
    }

    // -- Flatten real/imag parts depending on variant --
    vector<float> real_flat, imag_flat;

    auto flatten = [&](const auto& container) {
        using T = decay_t<decltype(container)>;

        if constexpr (is_same_v<T, cfloat>) {
            real_flat.push_back(container.real());
            if (is_complex) imag_flat.push_back(container.imag());
        }
        else if constexpr (is_same_v<T, vector<cfloat>>) {
            for (auto& v : container) {
                real_flat.push_back(v.real());
                if (is_complex) imag_flat.push_back(v.imag());
            }
        }
        else if constexpr (is_same_v<T, vector<vector<cfloat>>>) {
            for (auto& row : container) {
                for (auto& v : row) {
                    real_flat.push_back(v.real());
                    if (is_complex) imag_flat.push_back(v.imag());
                }
            }
        }
    };
    std::visit(flatten, data);

    // -- Write datasets under /values --
    hsize_t dims[1] = { real_flat.size() };
    DataSpace space(1, dims);

    file.createGroup("/values");
    // Create /values group if missing
    Group values_group = [&]() {
            return file.openGroup("/values");
    }();

    // Real part
    DataSet ds_real = values_group.createDataSet("real", PredType::NATIVE_FLOAT, space);
    ds_real.write(real_flat.data(), PredType::NATIVE_FLOAT);

    // Imag part (if complex)
    if (is_complex) {
        DataSet ds_imag = values_group.createDataSet("imag", PredType::NATIVE_FLOAT, space);
        ds_imag.write(imag_flat.data(), PredType::NATIVE_FLOAT);
    }
}
