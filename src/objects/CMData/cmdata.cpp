#include "H5Cpp.h"
#include "cmdata.hpp"
#include <algorithm>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>

using namespace H5;

CMData::CMData() {
    points = vector<Vec>();
    w_points = vector<float>();
    values = vector<complex<Vec>>();
    dimension = 0;
    is_complex = false;
    is_vector = false;
    with_w = false;
    with_n = false;
    n_inds = vector<int>();

    filled = false;
}

CMData::~CMData() {
    // Explicitly clear vectors (optional; only necessary if you're releasing large memory early)
    points.clear();
    w_points.clear();
    values.clear();
    n_inds.clear();

    // Reset flags (optional; mostly useful in debugging)
    filled = false;
}

CMData::CMData(vector<Vec> &points, vector<complex<Vec>> &values, int dimension,
               bool with_w, bool with_n, bool is_complex, bool is_vector) {
    this->points = points;
    this->values = values;
    this->dimension = dimension;
    this->is_complex = is_complex;
    this->is_vector = is_vector;
    this->with_w = with_w;
    this->with_n = with_n;
    this->n_inds = vector<int>();

    if (with_n)
        n_inds.push_back(0);
    if (with_w or with_n) {
        for (int i = 0; i < points.size(); i++) {
            if (with_w and find(w_points.begin(), w_points.end(),
                                points[i](dimension)) == w_points.end())
                w_points.push_back(points[i](dimension));
            if (with_n and i < points.size() - 1 and
                points[i + 1].n - points[i].n != 0)
                n_inds.push_back(i + 1);
        }
    }

    filled = true;
}

CMData::CMData(std::string filename) {
    using std::filesystem::exists;
    using std::filesystem::path;

    if (!exists(filename)) {
        throw std::runtime_error("File not found: " + filename);
    }

    std::string ext = path(filename).extension().string();

    if (ext == ".h5" || ext == ".hdf5") {
        load_hdf5(filename);
    } else if (ext == ".dat" || ext == ".txt") {
        // Assuming you already have load_dat defined somewhere
        *this = load(filename);  // if `load()` returns CMData
    } else {
        throw std::runtime_error("Unsupported file format: " + ext);
    }
}

CMData load(string filename) {
    ifstream file(filename);
    if (!file) {
        throw runtime_error("Failed to open file: " + filename);
    }

    string line;
    int dimension = 0;
    bool w_col = false, n_col = false, is_complex = false, is_vector = false;

    vector<Vec> points;
    vector<complex<Vec>> values;

    // Parse header
    if (!getline(file, line)) {
        throw runtime_error("Empty file or read error: " + filename);
    }

    istringstream header_stream(line);
    string word;
    while (header_stream >> word) {
        if (word == "x" || word == "kx" || word == "y" || word == "ky" || word == "z" || word == "kz")
            ++dimension;
        else if (word == "w") w_col = true;
        else if (word == "n") n_col = true;
        else if (word == "Re(f)") is_complex = true;
        else if (word == "fx") is_vector = true;
        else if (word == "Re(fx)") { is_vector = true; is_complex = true; }
    }

    const int coord_count = dimension + (w_col ? 1 : 0);
    float temp;

    // Read data lines
    while (getline(file, line)) {
        istringstream iss(line);
        Vec point;
        for (int i = 0; i < coord_count; ++i) {
            iss >> temp;
            point(i) = temp;
        }

        if (n_col) {
            iss >> temp;
            point.n = static_cast<int>(temp);
        }
        point.dimension = dimension;

        points.push_back(point);

        // Parse value(s)
        iss >> temp;
        Vec f_val(temp), fi_val(0.0f);

        if (is_complex) {
            iss >> temp;
            fi_val = Vec(temp);
        }

        if (is_vector) {
            f_val.dimension = dimension;
            fi_val.dimension = dimension;
            for (int i = 1; i < dimension; ++i) {
                iss >> temp;
                f_val(i) = temp;
                if (is_complex) {
                    iss >> temp;
                    fi_val(i) = temp;
                }
            }
        }

        values.emplace_back(f_val, fi_val);
    }

    return CMData(points, values, dimension, w_col, n_col, is_complex, is_vector);
}

string get_header(int dimension, bool with_w, bool with_n, bool is_complex,
                  bool is_vector) {
    // Determine fheader
    string fheader = " f    ";
    if (is_complex) {
        fheader = "   Re(f)         Im(f)";
    }
    if (is_vector && !is_complex) {
        fheader = "   fx    ";
    }
    if (is_vector && is_complex) {
        fheader = "    Re(fx)      Im(fx) ";
    }

    // Build header
    string header = "    x        ";
    if (dimension == 0) {
        header = "    w         ";
    }
    if (dimension > 1) {
        header += "    y         ";
        if (is_vector && !is_complex) {
            fheader += "   fy        ";
        } else if (is_vector && is_complex) {
            fheader += "    Re(fy)    Im(fy) ";
        }
    }
    if (dimension > 2) {
        header += "    z        ";
        if (is_vector && !is_complex) {
            fheader += "   fz        ";
        }
        if (is_vector && is_complex) {
            fheader += "     Re(fz)    Im(fz) ";
        }
    }
    if (with_w and dimension != 0) {
        header += "   w        ";
    }
    if (with_n) {
        header += "n        ";
    }

    header += fheader;
    return header;
}

void save_to_file(string filename, vector<Vec> &points,
                  vector<complex<Vec>> &values, int dimension, bool with_w,
                  bool with_n, bool is_complex, bool is_vector) {
    ofstream file(filename);
    string header =
        get_header(dimension, with_w, with_n, is_complex, is_vector);
    file << fixed << setprecision(6);
    file << header << endl;
    for (int i = 0; i < points.size(); i++) {
        Vec point = points[i];
        complex<Vec> value = values[i];
        if (dimension > 0) {
            if (point(0) >= 0)
                file << " ";
            file << point(0) << "    ";
            // if (point(0) < 0) file << " ";
        }

        if (dimension > 1) {
            if (point(1) >= 0)
                file << " ";
            file << point(1) << "    ";
            // if (point(1) < 0)
            //   file << " ";
        }

        if (dimension > 2) {
            if (point(2) >= 0)
                file << " ";
            file << point(2) << "    ";
            // if (point(2) < 0)
            //   file << " ";
        }

        if (with_w) {
            if (point(dimension) >= 0)
                file << " ";
            file << point(dimension) << "    ";
            // if (point(dimension) < 0)
            //   file << " ";
        }
        //if (point(dimension) >= 0)
        //    file << " ";
        if (with_n)
            file << point.n << "     ";

        if (value.real()(0) >= 0)
            file << " ";

        if (is_complex) {
            file << value.real()(0) << "     ";
            if (value.imag()(0) >= 0)
                file << " ";
            file << value.imag()(0) << "      ";
            if (is_vector and dimension > 1) {
                file << value.real()(1) << "      " << value.imag()(1)
                     << "      ";
                if (dimension > 2) {
                    file << value.real()(2) << "      " << value.imag()(2)
                         << "      ";
                }
            }
        } else {
            file << value.real()(0) << "      ";
            if (is_vector and dimension > 1) {
                file << value.real()(1) << "      ";
                if (dimension > 2) {
                    file << value.real()(2) << "      ";
                }
            }
        }
        file << endl;
    }
}

void save(CMData &data, string filename) {
    using std::filesystem::path;
    std::string ext = path(filename).extension().string();

    if (ext == ".h5" || ext == ".hdf5") {
        data.save_hdf5(filename);
    } else if (ext == ".dat" || ext == ".txt") {
        // Assuming you already have load_dat defined somewhere
        save_to_file(filename, data.points, data.values, data.dimension, data.with_w, data.with_n, data.is_complex, data.is_vector);
    } else {
        throw std::runtime_error("Unsupported file format: " + ext + "\nNot saving\n");
    }
}

void CMData::save_hdf5(const std::string& filename) const {
    H5File file(filename, H5F_ACC_TRUNC);

    // === Attributes ===
    file.createAttribute("dimension", PredType::NATIVE_INT, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_INT, &dimension);
    file.createAttribute("is_complex", PredType::NATIVE_HBOOL, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_HBOOL, &is_complex);
    file.createAttribute("is_vector",  PredType::NATIVE_HBOOL, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_HBOOL, &is_vector);
    file.createAttribute("with_w",     PredType::NATIVE_HBOOL, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_HBOOL, &with_w);
    file.createAttribute("with_n",     PredType::NATIVE_HBOOL, DataSpace(H5S_SCALAR)).write(PredType::NATIVE_HBOOL, &with_n);

    // === Points ===
    hsize_t n_points = points.size();

    hsize_t point_dims[2] = { n_points, static_cast<hsize_t>(dimension) };
    DataSpace point_space(2, point_dims);
    std::vector<float> point_data(n_points * dimension);
    for (size_t i = 0; i < n_points; ++i) {
        const Vec& v = points[i];
        float* dst = &point_data[i * dimension];
        if (dimension >= 1) dst[0] = v.x;
        if (dimension >= 2) dst[1] = v.y;
        if (dimension >= 3) dst[2] = v.z;
        if (dimension >= 4) dst[3] = v.w;
    }
    file.createDataSet("points", PredType::NATIVE_FLOAT, point_space).write(point_data.data(), PredType::NATIVE_FLOAT);

    // === Optional Weights ===
    if (with_w && !w_points.empty()) {
        hsize_t nw_points = w_points.size();
        hsize_t w_dims[1] = { nw_points };
        DataSpace w_space(1, w_dims);
        file.createDataSet("w_points", PredType::NATIVE_FLOAT, w_space).write(w_points.data(), PredType::NATIVE_FLOAT);
    }

    // === Optional N Indices ===
    if (with_n && !n_inds.empty()) {
        hsize_t n_dims[1] = { n_inds.size() };
        DataSpace n_space(1, n_dims);
        file.createDataSet("n_inds", PredType::NATIVE_INT, n_space).write(n_inds.data(), PredType::NATIVE_INT);
    }

    // === Values (complex<Vec>) ===
    int fdim = dimension;
    if (fdim == 0 and with_w) {
        fdim = 1;
    }
    hsize_t val_dims[3] = { n_points, 2, static_cast<hsize_t>(fdim) };
    DataSpace val_space(3, val_dims);
    std::vector<float> val_data(n_points * 2 * fdim);
    for (size_t i = 0; i < n_points; ++i) {
        const Vec& r = values[i].real();
        const Vec& im = values[i].imag();
        for (int d = 0; d < fdim; ++d) {
            val_data[i * 2 * fdim + 0 * fdim + d] = (&r.x)[d];
            val_data[i * 2 * fdim + 1 * fdim + d] = (&im.x)[d];
        }
    }
    file.createDataSet("values", PredType::NATIVE_FLOAT, val_space).write(val_data.data(), PredType::NATIVE_FLOAT);
}

void CMData::load_hdf5(const std::string& filename) {
    H5File file(filename, H5F_ACC_RDONLY);

    // === Read Attributes ===
    file.openAttribute("dimension").read(PredType::NATIVE_INT, &dimension);
    file.openAttribute("is_complex").read(PredType::NATIVE_HBOOL, &is_complex);
    file.openAttribute("is_vector").read(PredType::NATIVE_HBOOL, &is_vector);
    file.openAttribute("with_w").read(PredType::NATIVE_HBOOL, &with_w);
    file.openAttribute("with_n").read(PredType::NATIVE_HBOOL, &with_n);

    // === Read Points ===
    DataSet points_ds = file.openDataSet("points");
    DataSpace points_space = points_ds.getSpace();
    hsize_t dims[2];
    points_space.getSimpleExtentDims(dims);
    size_t n_points = dims[0];

    std::vector<float> point_data(n_points * dimension);
    points_ds.read(point_data.data(), PredType::NATIVE_FLOAT);
    points.resize(n_points);
    for (size_t i = 0; i < n_points; ++i) {
        Vec v = {};
        if (dimension >= 1) v.x = point_data[i * dimension + 0];
        if (dimension >= 2) v.y = point_data[i * dimension + 1];
        if (dimension >= 3) v.z = point_data[i * dimension + 2];
        if (dimension >= 4) v.w = point_data[i * dimension + 3];
        points[i] = v;
    }

    // === Optional Weights ===
    w_points.clear();
    if (with_w) {
        try {
            DataSet w_ds = file.openDataSet("w_points");
            hsize_t w_dims[1];
            w_ds.getSpace().getSimpleExtentDims(w_dims);
            w_points.resize(w_dims[0]);
            w_ds.read(w_points.data(), PredType::NATIVE_FLOAT);
        } catch (...) {
            printf("Skipping. Error with with_w\n");
            // Dataset missing: silently skip
        }
    }

    // === Optional N Indices ===
    n_inds.clear();
    if (with_n) {
        try {
            DataSet n_ds = file.openDataSet("n_inds");
            hsize_t n_dims[1];
            n_ds.getSpace().getSimpleExtentDims(n_dims);
            n_inds.resize(n_dims[0]);
            n_ds.read(n_inds.data(), PredType::NATIVE_INT);
        } catch (...) {
            printf("Skipping. Error with with_n\n");
            // Dataset missing: silently skip
        }
    }

    // === Read Values ===
    DataSet val_ds = file.openDataSet("values");
    DataSpace val_space = val_ds.getSpace();
    hsize_t val_dims[3];
    val_space.getSimpleExtentDims(val_dims);
    std::vector<float> val_data(val_dims[0] * val_dims[1] * val_dims[2]);
    val_ds.read(val_data.data(), PredType::NATIVE_FLOAT);

    values.resize(n_points);
    int fdim = dimension;
    if (fdim == 0 and with_w)
        fdim = 1;
    for (size_t i = 0; i < n_points; ++i) {
        Vec r = {}, im = {};
        for (int d = 0; d < fdim; ++d) {
            (&r.x)[d] = val_data[i * 2 * fdim + 0 * fdim + d];
            (&im.x)[d] = val_data[i * 2 * fdim + 1 * fdim + d];
        }
        values[i] = std::complex<Vec>(r, im);
    }

    filled = true;
}
