// base_field.cpp
#include "base_data.hpp"
#include <H5Cpp.h>
#include <stdexcept>
#include <iostream>
#include <cmath>

using namespace H5;

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
    // Read bools as ints to avoid size mismatch issues
    int temp_is_complex, temp_is_vector, temp_is_matrix, temp_with_k, temp_with_w, temp_as_mesh;
    file.openDataSet("/is_complex").read(&temp_is_complex, PredType::NATIVE_INT);
    file.openDataSet("/is_vector").read(&temp_is_vector, PredType::NATIVE_INT);
    file.openDataSet("/is_matrix").read(&temp_is_matrix, PredType::NATIVE_INT);
    file.openDataSet("/with_k").read(&temp_with_k, PredType::NATIVE_INT);
    file.openDataSet("/with_w").read(&temp_with_w, PredType::NATIVE_INT);

    field.is_complex = temp_is_complex;
    field.is_vector = temp_is_vector;
    field.is_matrix = temp_is_matrix;
    field.with_k = temp_with_k;
    field.with_w = temp_with_w;

    // Read inds (tensor indices)
    try {
        DataSet ds_inds = file.openDataSet("/inds");
        DataSpace space_inds = ds_inds.getSpace();
        hsize_t inds_size;
        space_inds.getSimpleExtentDims(&inds_size);
        field.inds.resize(inds_size);
        ds_inds.read(field.inds.data(), PredType::NATIVE_INT);
        space_inds.close();
        ds_inds.close();
    } catch (...) {
        // Legacy format: use empty inds for scalar
        field.inds = {};
    }

    file.openDataSet("/dimension").read(&field.dimension, PredType::NATIVE_INT);

    file.openDataSet("/as_mesh").read(&temp_as_mesh, PredType::NATIVE_INT);
    field.as_mesh = temp_as_mesh;

    // -- Domain (optional) --
    if (field.with_k) {
        DataSet ds_domain = file.openDataSet("/domain");
        DataSpace space_domain = ds_domain.getSpace();
        int rank = space_domain.getSimpleExtentNdims();
        if (rank == 2) {
            hsize_t dims[2];
            space_domain.getSimpleExtentDims(dims);
            // Read into flat buffer first, then populate 2D structure
            std::vector<float> domain_flat(dims[0] * dims[1]);
            ds_domain.read(domain_flat.data(), PredType::NATIVE_FLOAT);
            field.domain.assign(dims[0], std::vector<float>(dims[1]));
            for (size_t i = 0; i < dims[0]; i++) {
                for (size_t j = 0; j < dims[1]; j++) {
                    field.domain[i][j] = domain_flat[i * dims[1] + j];
                }
            }
        }
        space_domain.close();
        ds_domain.close();
    }

    // -- Mesh (optional) --
    if (field.as_mesh) {
        DataSet ds = file.openDataSet("/mesh");
        DataSpace space = ds.getSpace();
        hsize_t dims[1];
        space.getSimpleExtentDims(dims);
        field.mesh.resize(dims[0]);
        ds.read(field.mesh.data(), PredType::NATIVE_INT);
        space.close();
        ds.close();
    }

    // -- Points (k-point data when as_mesh = false) --
    if (!field.as_mesh && field.with_k) {
        try {
            DataSet ds = file.openDataSet("/points");
            DataSpace space = ds.getSpace();
            int rank = space.getSimpleExtentNdims();
            if (rank == 2) {
                hsize_t dims[2];
                space.getSimpleExtentDims(dims);
                // Read into flat buffer first, then populate 2D structure
                std::vector<float> points_flat(dims[0] * dims[1]);
                ds.read(points_flat.data(), PredType::NATIVE_FLOAT);
                field.points.assign(dims[0], std::vector<float>(dims[1]));
                for (size_t i = 0; i < dims[0]; i++) {
                    for (size_t j = 0; j < dims[1]; j++) {
                        field.points[i][j] = points_flat[i * dims[1] + j];
                    }
                }
            }
            space.close();
            ds.close();
        } catch (...) {
            // ignore if missing
        }
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
            space.close();
            ds.close();
        } catch (...) {
            // ignore if missing
        }
    }

    // -- Compute sizes --
    int total_indices = field.total_index_size();  // FIX: Use new method instead of buggy pow
    int nk = field.nk();
    int nw = field.nw();
    int vec_len = field.vec_len();
    int total_elements = total_indices * nk * nw * vec_len;

    // --- Real part ---
    H5::DataSet real_ds = file.openDataSet("/values/real");
    H5::DataSpace real_space = real_ds.getSpace();
    hsize_t dims[1];
    real_space.getSimpleExtentDims(dims);

    std::vector<float> real_flat(dims[0]);
    real_ds.read(real_flat.data(), H5::PredType::NATIVE_FLOAT);
    real_space.close();
    real_ds.close();

    // --- Imag part (if complex) ---
    std::vector<float> imag_flat;
    if (field.is_complex) {
        H5::DataSet imag_ds = file.openDataSet("/values/imag");
        H5::DataSpace imag_space = imag_ds.getSpace();
        hsize_t dims_imag[1];
        imag_space.getSimpleExtentDims(dims_imag);

        imag_flat.resize(dims_imag[0]);
        imag_ds.read(imag_flat.data(), H5::PredType::NATIVE_FLOAT);
        imag_space.close();
        imag_ds.close();
    }

    // -- Populate variant based on tensor rank --
    int rank = field.rank();
    int total_tensors = nk * nw;

    if (rank == 4) {
        // 4D tensor: [nk*nw][d0][d1][d2][d3] with potentially non-uniform dimensions
        int d0 = field.inds[0], d1 = field.inds[1], d2 = field.inds[2], d3 = field.inds[3];
        std::vector<std::vector<std::vector<std::vector<std::vector<cfloat>>>>> tensors(total_tensors);

        int idx = 0;
        for (int t = 0; t < total_tensors; ++t) {
            tensors[t].resize(d0);
            for (int i = 0; i < d0; ++i) {
                tensors[t][i].resize(d1);
                for (int j = 0; j < d1; ++j) {
                    tensors[t][i][j].resize(d2);
                    for (int k = 0; k < d2; ++k) {
                        tensors[t][i][j][k].resize(d3);
                        for (int l = 0; l < d3; ++l) {
                            tensors[t][i][j][k][l] = cfloat(real_flat[idx], field.is_complex ? imag_flat[idx] : 0.0f);
                            idx++;
                        }
                    }
                }
            }
        }
        field.data = tensors;
    } else if (rank == 3) {
        // 3D tensor: [nk*nw][d0][d1][d2]
        int d0 = field.inds[0], d1 = field.inds[1], d2 = field.inds[2];
        std::vector<std::vector<std::vector<std::vector<cfloat>>>> tensors(total_tensors);

        int idx = 0;
        for (int t = 0; t < total_tensors; ++t) {
            tensors[t].resize(d0);
            for (int i = 0; i < d0; ++i) {
                tensors[t][i].resize(d1);
                for (int j = 0; j < d1; ++j) {
                    tensors[t][i][j].resize(d2);
                    for (int k = 0; k < d2; ++k) {
                        tensors[t][i][j][k] = cfloat(real_flat[idx], field.is_complex ? imag_flat[idx] : 0.0f);
                        idx++;
                    }
                }
            }
        }
        field.data = tensors;
    } else if (rank == 2) {
        // Matrix: [nk*nw][d0][d1] (supports non-square matrices!)
        int d0 = field.inds[0], d1 = field.inds[1];
        std::vector<std::vector<std::vector<cfloat>>> matrices(total_tensors);

        int idx = 0;
        for (int t = 0; t < total_tensors; ++t) {
            matrices[t].resize(d0);
            for (int i = 0; i < d0; ++i) {
                matrices[t][i].resize(d1);
                for (int j = 0; j < d1; ++j) {
                    matrices[t][i][j] = cfloat(real_flat[idx], field.is_complex ? imag_flat[idx] : 0.0f);
                    idx++;
                }
            }
        }
        field.data = matrices;
    } else if (rank == 1) {
        // Vector: [nk*nw][d0]
        int d0 = field.inds[0];
        std::vector<std::vector<cfloat>> vectors(total_tensors);

        int idx = 0;
        for (int t = 0; t < total_tensors; ++t) {
            vectors[t].resize(d0);
            for (int i = 0; i < d0; ++i) {
                vectors[t][i] = cfloat(real_flat[idx], field.is_complex ? imag_flat[idx] : 0.0f);
                idx++;
            }
        }
        field.data = vectors;
    } else {
        // Scalar: rank = 0, inds = {}
        std::vector<cfloat> scalars(total_tensors);
        for (int i = 0; i < total_tensors; ++i) {
            scalars[i] = cfloat(real_flat[i], field.is_complex ? imag_flat[i] : 0.0f);
        }
        field.data = scalars;
    }

    // Explicitly close the file to release locks immediately
    file.close();

    return field;
}

// Load with specified ordering (k-w or w-k)
BaseData load_data_from_hdf5(const std::string& filename, const std::string& ordering) {
    BaseData field = load_data_from_hdf5(filename);

    // Check if we need to reorder (stored as k-w by default)
    if (ordering == "w-k" && field.with_k && field.with_w) {
        int nk = field.nk();
        int nw = field.nw();

        // Reorder based on tensor rank
        int rank = field.rank();

        if (rank == 4) {
            // 4D tensor data
            auto& tensors = field.get<std::vector<std::vector<std::vector<std::vector<std::vector<cfloat>>>>>>();
            std::vector<std::vector<std::vector<std::vector<std::vector<cfloat>>>>> reordered(nk * nw);

            for (int k = 0; k < nk; ++k) {
                for (int w = 0; w < nw; ++w) {
                    reordered[w * nk + k] = tensors[k * nw + w];
                }
            }
            field.data = reordered;
        } else if (rank == 3) {
            // 3D tensor data
            auto& tensors = field.get<std::vector<std::vector<std::vector<std::vector<cfloat>>>>>();
            std::vector<std::vector<std::vector<std::vector<cfloat>>>> reordered(nk * nw);

            for (int k = 0; k < nk; ++k) {
                for (int w = 0; w < nw; ++w) {
                    reordered[w * nk + k] = tensors[k * nw + w];
                }
            }
            field.data = reordered;
        } else if (rank == 2) {
            // Matrix data
            auto& matrices = field.get<std::vector<std::vector<std::vector<cfloat>>>>();
            std::vector<std::vector<std::vector<cfloat>>> reordered(nk * nw);

            for (int k = 0; k < nk; ++k) {
                for (int w = 0; w < nw; ++w) {
                    reordered[w * nk + k] = matrices[k * nw + w];
                }
            }
            field.data = reordered;
        } else if (rank == 1) {
            // Vector data
            auto& vectors = field.get<std::vector<std::vector<cfloat>>>();
            std::vector<std::vector<cfloat>> reordered(nk * nw);

            for (int k = 0; k < nk; ++k) {
                for (int w = 0; w < nw; ++w) {
                    reordered[w * nk + k] = vectors[k * nw + w];
                }
            }
            field.data = reordered;
        } else {
            // Scalar data
            auto& scalars = field.get<std::vector<cfloat>>();
            std::vector<cfloat> reordered(nk * nw);

            for (int k = 0; k < nk; ++k) {
                for (int w = 0; w < nw; ++w) {
                    reordered[w * nk + k] = scalars[k * nw + w];
                }
            }
            field.data = reordered;
        }
    }

    return field;
}

void save_data_to_hdf5(BaseData& field, const std::string& filename) {
    save_data_to_hdf5(filename,
                       field.is_complex,
                       field.is_vector,
                       field.is_matrix,
                       field.with_k,
                       field.with_w,
                       field.as_mesh,
                       field.inds,  // CHANGED: use inds instead of n_indices/dim_indices
                       field.mesh,
                       field.domain,
                       field.dimension,
                       field.w_points,
                       field.points,
                       field.data);
}

// Save with specified ordering (k-w or w-k)
void save_data_to_hdf5(BaseData& field, const std::string& filename, const std::string& ordering) {
    // If ordering is k-w or default, just save normally
    if (ordering == "k-w" || ordering.empty()) {
        save_data_to_hdf5(field, filename);
        return;
    }

    // If ordering is w-k, we need to reorder the data before saving
    if (ordering == "w-k" && field.with_k && field.with_w) {
        BaseData reordered = field;  // Copy
        int nk = field.nk();
        int nw = field.nw();

        // Reorder based on tensor rank (same as in load function)
        int rank = field.rank();

        if (rank == 4) {
            auto& tensors = field.get<std::vector<std::vector<std::vector<std::vector<std::vector<cfloat>>>>>>();
            std::vector<std::vector<std::vector<std::vector<std::vector<cfloat>>>>> reordered_data(nk * nw);
            for (int k = 0; k < nk; ++k) {
                for (int w = 0; w < nw; ++w) {
                    reordered_data[w * nk + k] = tensors[k * nw + w];
                }
            }
            reordered.data = reordered_data;
        } else if (rank == 3) {
            auto& tensors = field.get<std::vector<std::vector<std::vector<std::vector<cfloat>>>>>();
            std::vector<std::vector<std::vector<std::vector<cfloat>>>> reordered_data(nk * nw);
            for (int k = 0; k < nk; ++k) {
                for (int w = 0; w < nw; ++w) {
                    reordered_data[w * nk + k] = tensors[k * nw + w];
                }
            }
            reordered.data = reordered_data;
        } else if (rank == 2) {
            auto& matrices = field.get<std::vector<std::vector<std::vector<cfloat>>>>();
            std::vector<std::vector<std::vector<cfloat>>> reordered_data(nk * nw);
            for (int k = 0; k < nk; ++k) {
                for (int w = 0; w < nw; ++w) {
                    reordered_data[w * nk + k] = matrices[k * nw + w];
                }
            }
            reordered.data = reordered_data;
        } else if (rank == 1) {
            auto& vectors = field.get<std::vector<std::vector<cfloat>>>();
            std::vector<std::vector<cfloat>> reordered_data(nk * nw);
            for (int k = 0; k < nk; ++k) {
                for (int w = 0; w < nw; ++w) {
                    reordered_data[w * nk + k] = vectors[k * nw + w];
                }
            }
            reordered.data = reordered_data;
        } else {
            auto& scalars = field.get<std::vector<cfloat>>();
            std::vector<cfloat> reordered_data(nk * nw);
            for (int k = 0; k < nk; ++k) {
                for (int w = 0; w < nw; ++w) {
                    reordered_data[w * nk + k] = scalars[k * nw + w];
                }
            }
            reordered.data = reordered_data;
        }

        save_data_to_hdf5(reordered, filename);
    } else {
        save_data_to_hdf5(field, filename);
    }
}

void save_data(string filename, BaseData::DataVariant& data, bool is_complex, vector<int> mesh, vector<vector<float>> domain, vector<float> w_points, const vector<int>& inds) {
    bool is_vector = false; // Will add vector support when it becomes relevant
    bool with_k = mesh.size() > 0;
    bool with_w = w_points.size() > 0;
    bool as_mesh = mesh.size() > 0; // Would be false if points were given or if mesh is empty

    // Determine dimension from mesh: count non-trivial dimensions (mesh[i] > 1)
    int dim = 0;
    for (size_t i = 0; i < mesh.size(); i++) {
        if (mesh[i] > 1) dim++;
    }
    if (dim == 0) dim = domain.empty() ? 3 : domain.size();  // Fallback to domain size

    bool is_matrix = inds.size() == 2;  // Matrix if rank = 2
    vector<vector<float>> points = {}; // Empty for this wrapper function
    save_data_to_hdf5(filename, is_complex, is_vector, is_matrix, with_k, with_w, as_mesh, inds, mesh, domain, dim, w_points, points, data);
}

void save_data_to_hdf5(const std::string& filename,
                        bool is_complex,
                        bool is_vector,
                        bool is_matrix,
                        bool with_k,
                        bool with_w,
                        bool as_mesh,
                        const std::vector<int>& inds,  // CHANGED: inds instead of n_indices/dim_indices
                        std::vector<int>& mesh,
                        std::vector<std::vector<float>>& domain,
                        int dimension,
                        std::vector<float>& w_points,
                        std::vector<std::vector<float>>& points,
                        const BaseData::DataVariant& data) {
    H5File file(filename, H5F_ACC_TRUNC);

    // -- Metadata scalars --
    auto write_scalar = [&](const std::string& name, int value) {
        DataSpace scalar_space(H5S_SCALAR);
        DataSet ds = file.createDataSet(name, PredType::NATIVE_INT, scalar_space);
        ds.write(&value, PredType::NATIVE_INT);
    };

    write_scalar("/is_complex", is_complex);
    write_scalar("/is_vector", is_vector);
    write_scalar("/is_matrix", is_matrix);
    write_scalar("/with_k", with_k);
    write_scalar("/with_w", with_w);
    write_scalar("/dimension", dimension);
    write_scalar("/as_mesh", as_mesh);

    // -- Write inds (tensor indices) as array --
    if (!inds.empty()) {
        hsize_t dims[1] = { inds.size() };
        DataSpace space(1, dims);
        DataSet ds = file.createDataSet("/inds", PredType::NATIVE_INT, space);
        ds.write(inds.data(), PredType::NATIVE_INT);
    } else {
        // Write empty inds for scalar
        hsize_t dims[1] = { 0 };
        DataSpace space(1, dims);
        DataSet ds = file.createDataSet("/inds", PredType::NATIVE_INT, space);
    }

    // -- Mesh --
    if (!mesh.empty()) {
        hsize_t dims[1] = { mesh.size() };
        DataSpace space(1, dims);
        DataSet ds = file.createDataSet("/mesh", PredType::NATIVE_INT, space);
        ds.write(mesh.data(), PredType::NATIVE_INT);
    }

    // -- Points (k-point data when as_mesh = false) --
    if (!points.empty()) {
        hsize_t dims[2] = { points.size(), points[0].size() };
        DataSpace space(2, dims);
        DataSet ds = file.createDataSet("/points", PredType::NATIVE_FLOAT, space);

        // Flatten 2D into 1D buffer
        std::vector<float> flat;
        flat.reserve(points.size() * points[0].size());
        for (auto const& point : points) {
            flat.insert(flat.end(), point.begin(), point.end());
        }
        ds.write(flat.data(), PredType::NATIVE_FLOAT);
    }

    // -- Domain --
    if (!domain.empty()) {
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

    // -- w_points --
    if (!w_points.empty()) {
        hsize_t dims[1] = { w_points.size() };
        DataSpace space(1, dims);
        DataSet ds = file.createDataSet("/w_points", PredType::NATIVE_FLOAT, space);
        ds.write(w_points.data(), PredType::NATIVE_FLOAT);
    }
    else {
        hsize_t dims[1] = { w_points.size() };
        DataSpace space(1, dims);
        DataSet ds = file.createDataSet("/w_points", PredType::NATIVE_FLOAT, space);
        ds.write(w_points.data(), PredType::NATIVE_FLOAT);
    }

    // -- Flatten real/imag parts depending on variant --
    std::vector<float> real_flat, imag_flat;

    auto flatten = [&](auto const& container) {
        using T = std::decay_t<decltype(container)>;
        if constexpr (std::is_same_v<T, std::vector<cfloat>>) {
            // Scalar
            for (auto& v : container) {
                real_flat.push_back(v.real());
                if (is_complex) imag_flat.push_back(v.imag());
            }
        } else if constexpr (std::is_same_v<T, std::vector<std::vector<cfloat>>>) {
            // Vector (rank=1)
            for (auto& vec : container) {
                for (auto& v : vec) {
                    real_flat.push_back(v.real());
                    if (is_complex) imag_flat.push_back(v.imag());
                }
            }
        } else if constexpr (std::is_same_v<T, std::vector<std::vector<std::vector<cfloat>>>>) {
            // Matrix (rank=2)
            for (auto& matrix : container) {
                for (auto& row : matrix) {
                    for (auto& v : row) {
                        real_flat.push_back(v.real());
                        if (is_complex) imag_flat.push_back(v.imag());
                    }
                }
            }
        } else if constexpr (std::is_same_v<T, std::vector<std::vector<std::vector<std::vector<cfloat>>>>>) {
            // 3D tensor (rank=3)
            for (auto& tensor : container) {
                for (auto& plane : tensor) {
                    for (auto& row : plane) {
                        for (auto& v : row) {
                            real_flat.push_back(v.real());
                            if (is_complex) imag_flat.push_back(v.imag());
                        }
                    }
                }
            }
        } else if constexpr (std::is_same_v<T, std::vector<std::vector<std::vector<std::vector<std::vector<cfloat>>>>>>) {
            // 4D tensor (rank=4)
            for (auto& tensor : container) {
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
        }
    };
    std::visit(flatten, data);

    // -- Write datasets under /values --
    hsize_t dims[1] = { real_flat.size() };
    DataSpace space(1, dims);

    // Ensure /values group exists
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
