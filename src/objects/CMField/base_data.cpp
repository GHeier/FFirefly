// base_field.cpp
#include "base_data.hpp"
#include <H5Cpp.h>
#include <stdexcept>
#include <iostream>
#include <cmath>

using namespace H5;

static map<string, string> read_attributes(Group& root, int& format_version);

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
    field.centered = true;  // Default for legacy files
    if (file.nameExists("/centered")) {
        int temp_centered;
        file.openDataSet("/centered").read(&temp_centered, PredType::NATIVE_INT);
        field.centered = temp_centered;
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

int infer_dimension(const vector<int>& mesh, const vector<vector<float>>& points) {
    if (mesh.empty())
        return points.empty() ? 0 : points[0].size();
    int dim = mesh.size();
    while (dim > 1 && mesh[dim - 1] == 1) dim--;
    return dim;
}

static bool domain_is_empty(const std::vector<std::vector<float>>& domain) {
    for (auto const& row : domain)
        if (!row.empty()) return false;
    return true;
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
    if (dims[0] == 0 || dims[1] == 0) {
        // Legacy files may store the default {{}} domain as a 1x0 dataset
        field.domain.clear();
        return;
    }
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
    {
        Group root = file.openGroup("/");
        field.attrs = read_attributes(root, field.format_version);
    }

    // -- Domain (optional) --
    if (field.with_k && field.as_mesh && file.nameExists("/domain")) {
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
// Flattening
// ============================================================================

namespace {

void append(const cfloat& v, vector<float>& re, vector<float>& im) {
    re.push_back(v.real());
    im.push_back(v.imag());
}
void append(float v, vector<float>& re, vector<float>& im) {
    re.push_back(v);
    im.push_back(0.0f);
}
template <typename T>
void append(const vector<T>& v, vector<float>& re, vector<float>& im) {
    for (auto& e : v) append(e, re, im);
}

bool variant_empty(const BaseData::DataVariant& d) {
    return std::visit([](auto& v) { return v.empty(); }, d);
}

} // namespace

// The structured variant wins when it holds data (it is what programmatic
// callers set); otherwise the flat arrays from HDF5 are used.
void BaseData::flat_values(vector<float>& re, vector<float>& im) const {
    const BaseData& f = *this;
    re.clear();
    im.clear();
    if (!variant_empty(f.data)) {
        std::visit([&](auto& d) { append(d, re, im); }, f.data);
    } else {
        re = f.real_values;
        im = f.imag_values;
        im.resize(re.size(), 0.0f);
    }
    if (!f.is_complex) im.clear();
}

void BaseData::flatten_variant() {
    if (!real_values.empty() || variant_empty(data)) return;
    flat_values(real_values, imag_values);
}

// ============================================================================
// Attributes
// ============================================================================

static void write_string_attribute(Group& root, const string& key, const string& value) {
    StrType type(PredType::C_S1, H5T_VARIABLE);
    if (root.attrExists(key)) root.removeAttr(key);
    Attribute attr = root.createAttribute(key, type, DataSpace(H5S_SCALAR));
    attr.write(type, value);
}

// Reads string, integer and float attributes; numbers are converted to text
static map<string, string> read_attributes(Group& root, int& format_version) {
    map<string, string> out;
    format_version = 1;
    int n = root.getNumAttrs();
    for (int i = 0; i < n; i++) {
        Attribute attr = root.openAttribute((unsigned)i);
        string name = attr.getName();
        H5T_class_t cls = attr.getTypeClass();
        if (name == "format_version" && cls == H5T_INTEGER) {
            attr.read(PredType::NATIVE_INT, &format_version);
        } else if (cls == H5T_STRING) {
            string value;
            attr.read(attr.getStrType(), value);
            out[name] = value;
        } else if (cls == H5T_INTEGER && attr.getSpace().getSimpleExtentNpoints() == 1) {
            long long value;
            attr.read(PredType::NATIVE_LLONG, &value);
            out[name] = std::to_string(value);
        } else if (cls == H5T_FLOAT && attr.getSpace().getSimpleExtentNpoints() == 1) {
            double value;
            attr.read(PredType::NATIVE_DOUBLE, &value);
            char buf[64];
            snprintf(buf, sizeof(buf), "%.9g", value);
            out[name] = buf;
        }
    }
    return out;
}

map<string, string> read_field_attributes(const std::string& filename) {
    H5File file(filename, H5F_ACC_RDONLY);
    Group root = file.openGroup("/");
    int version;
    return read_attributes(root, version);
}

void write_field_attributes(const std::string& filename, const map<string, string>& attrs) {
    H5File file(filename, H5F_ACC_RDWR);
    Group root = file.openGroup("/");
    for (auto& [key, value] : attrs) {
        if (key == "format_version")
            throw std::runtime_error("format_version is set by the library and cannot be overridden");
        write_string_attribute(root, key, value);
    }
}

// ============================================================================
// Saving
// ============================================================================

void save_data_to_hdf5(BaseData& field, const std::string& filename) {
    vector<float> real_flat, imag_flat;
    field.flat_values(real_flat, imag_flat);

    H5File file(filename, H5F_ACC_TRUNC);
    write_metadata(file, field.is_complex, field.is_vector, field.with_k, field.with_w,
                   field.as_mesh, field.centered, field.dimension, field.inds, field.mesh);

    if (!field.points.empty()) write_points(file, field.points);
    if (!domain_is_empty(field.domain)) write_domain(file, field.domain);
    if (!field.w_points.empty()) write_w_points(file, field.w_points);

    hsize_t dims[1] = { real_flat.size() };
    DataSpace space(1, dims);
    Group values_group = file.createGroup("/values");
    DataSet ds_real = values_group.createDataSet("real", PredType::NATIVE_FLOAT, space);
    ds_real.write(real_flat.data(), PredType::NATIVE_FLOAT);
    if (field.is_complex) {
        DataSet ds_imag = values_group.createDataSet("imag", PredType::NATIVE_FLOAT, space);
        ds_imag.write(imag_flat.data(), PredType::NATIVE_FLOAT);
    }

    Group root = file.openGroup("/");
    int version = BaseData::FORMAT_VERSION;
    Attribute va = root.createAttribute("format_version", PredType::NATIVE_INT, DataSpace(H5S_SCALAR));
    va.write(PredType::NATIVE_INT, &version);
    for (auto& [key, value] : field.attrs)
        if (key != "format_version") write_string_attribute(root, key, value);
}

// ============================================================================
// Convenience save_data wrappers
// ============================================================================

template <typename T>
static void save_nested(const string& filename, const T& data, bool is_complex,
                        const vector<int>& inds, const vector<int>& mesh,
                        const vector<vector<float>>& domain, const vector<float>& w_points,
                        const vector<vector<float>>& points, bool centered) {
    BaseData f;
    f.is_complex = is_complex;
    f.with_k = points.size() > 0 || mesh.size() > 0;
    f.with_w = w_points.size() > 0;
    f.as_mesh = mesh.size() > 0;
    f.centered = centered;
    f.inds = inds;
    f.mesh = mesh;
    f.domain = domain;
    f.w_points = w_points;
    f.points = points;
    f.dimension = infer_dimension(mesh, points);
    append(data, f.real_values, f.imag_values);
    if (!is_complex) f.imag_values.clear();
    save_data_to_hdf5(f, filename);
}

#define FIREFLY_SAVE_DATA(TYPE, IS_COMPLEX)                                                     \
    void save_data(string filename, TYPE& data, vector<int> inds, vector<int> mesh,             \
                   vector<vector<float>> domain, vector<float> w_points,                         \
                   const vector<vector<float>>& points, bool centered) {                         \
        save_nested(filename, data, IS_COMPLEX, inds, mesh, domain, w_points, points, centered); \
    }

FIREFLY_SAVE_DATA(vector<vector<vector<vector<cfloat>>>>, true)
FIREFLY_SAVE_DATA(vector<vector<vector<cfloat>>>, true)
FIREFLY_SAVE_DATA(vector<vector<cfloat>>, true)
FIREFLY_SAVE_DATA(vector<cfloat>, true)
FIREFLY_SAVE_DATA(vector<vector<vector<vector<float>>>>, false)
FIREFLY_SAVE_DATA(vector<vector<vector<float>>>, false)
FIREFLY_SAVE_DATA(vector<vector<float>>, false)
FIREFLY_SAVE_DATA(vector<float>, false)
#undef FIREFLY_SAVE_DATA

void save_data(string filename, vector<float>& data, vector<float> w_points, bool centered) {
    save_nested(filename, data, false, {}, {}, {}, w_points, {}, centered);
}
