#include "fields.hpp"
#include "field.hpp"
#include "../vec.hpp"
#include <openblas/lapacke.h>

// Field_C implementation
Field_C::Field_C()
    : cmf({}, true, false, false, {}, {}, {}) {}

Field_C::Field_C(const BaseData::DataVariant& data,
                 const vector<int>& mesh,
                 const vector<vector<float>>& domain,
                 const vector<float>& w_points)
    : cmf(data, true, false, false, mesh, domain, w_points) {}

Field_C::Field_C(FieldImpl f) : cmf(f) {}

Field_C::Field_C(const string& filename) : cmf(filename) {}

complex<float> Field_C::operator()(float w) {
    auto result = cmf(w);
    if (auto* c = std::get_if<cfloat>(&result)) {
        return *c;
    }
    // Shouldn't reach here for complex scalar field
    return complex<float>(0, 0);
}

complex<float> Field_C::operator()(Vec point, float w) {
    auto result = cmf(point, w);
    if (auto* c = std::get_if<cfloat>(&result)) {
        return *c;
    }
    // Shouldn't reach here for complex scalar field
    return complex<float>(0, 0);
}

Field_C& Field_C::operator=(const Field_C& other) {
    if (this != &other) {
        cmf = other.cmf;
    }
    return *this;
}

void Field_C::save(const string& filename) {
    cmf.save(filename);
}

vector<complex<float>> Field_C::operator()(const vector<Vec>& points, float w) {
    vector<complex<float>> results(points.size());
    for (size_t i = 0; i < points.size(); i++) {
        results[i] = (*this)(points[i], w);
    }
    return results;
}

vector<complex<float>> Field_C::operator()(const vector<float>& w_points) {
    vector<complex<float>> results(w_points.size());
    for (size_t i = 0; i < w_points.size(); i++) {
        results[i] = (*this)(w_points[i]);
    }
    return results;
}

BaseData* Field_C::get_data() {
    return &cmf.data;
}


// Field_R implementation
Field_R::Field_R()
    : cmf({}, false, false, false, {}, {}, {}) {}

Field_R::Field_R(const BaseData::DataVariant& data,
                 const vector<int>& mesh,
                 const vector<vector<float>>& domain,
                 const vector<float>& w_points)
    : cmf(data, false, false, false, mesh, domain, w_points) {}

Field_R::Field_R(FieldImpl f) : cmf(f) {}

Field_R::Field_R(const string& filename) : cmf(filename) {}

float Field_R::operator()(float w) {
    auto result = cmf(w);
    if (auto* r = std::get_if<float>(&result)) {
        return *r;
    }
    // Shouldn't reach here for real scalar field
    return 0.0f;
}

float Field_R::operator()(Vec point, float w) {
    auto result = cmf(point, w);
    if (auto* r = std::get_if<float>(&result)) {
        return *r;
    }
    // Shouldn't reach here for real scalar field
    return 0.0f;
}

Field_R& Field_R::operator=(const Field_R& other) {
    if (this != &other) {
        cmf = other.cmf;
    }
    return *this;
}

void Field_R::save(const string& filename) {
    cmf.save(filename);
}

vector<float> Field_R::operator()(const vector<Vec>& points, float w) {
    vector<float> results(points.size());
    for (size_t i = 0; i < points.size(); i++) {
        results[i] = (*this)(points[i], w);
    }
    return results;
}

vector<float> Field_R::operator()(const vector<float>& w_points) {
    vector<float> results(w_points.size());
    for (size_t i = 0; i < w_points.size(); i++) {
        results[i] = (*this)(w_points[i]);
    }
    return results;
}

BaseData* Field_R::get_data() {
    return &cmf.data;
}

// Field_CM implementation (Complex Matrix)
Field_CM::Field_CM()
    : cmf(vector<vector<vector<cfloat>>>(), true, false, true, {}, {}, {}, {3, 3}) {}

Field_CM::Field_CM(const BaseData::DataVariant& data,
                   const vector<int>& inds,
                   const vector<int>& mesh,
                   const vector<vector<float>>& domain,
                   const vector<float>& w_points)
    : cmf(data, true, false, true, mesh, domain, w_points, inds) {}

Field_CM::Field_CM(FieldImpl f) : cmf(f) {}

Field_CM::Field_CM(const string& filename) : cmf(filename) {}

Field_CM& Field_CM::operator=(const Field_CM& other) {
    if (this != &other) {
        cmf = other.cmf;
    }
    return *this;
}

void Field_CM::save(const string& filename) {
    cmf.save(filename);
}

vector<vector<cfloat>> Field_CM::operator()(Vec point, float w) {
    auto result = cmf.get_array(point, w);

    // Try 1D vector (rank=1) - wrap as column matrix
    if (auto* vec = std::get_if<vector<cfloat>>(&result)) {
        vector<vector<cfloat>> mat(vec->size(), vector<cfloat>(1));
        for (size_t i = 0; i < vec->size(); i++) {
            mat[i][0] = (*vec)[i];
        }
        return mat;
    }

    // Try 2D matrix (rank=2)
    if (auto* mat = std::get_if<vector<vector<cfloat>>>(&result)) {
        return *mat;
    }

    // Try 3D tensor (rank=3) - flatten to matrix by combining first two indices
    if (auto* ten3 = std::get_if<vector<vector<vector<cfloat>>>>(&result)) {
        if (ten3->empty()) return vector<vector<cfloat>>();
        int d1 = ten3->size();
        int d2 = (*ten3)[0].size();
        int d3 = (*ten3)[0][0].size();
        // Flatten [d1][d2][d3] to [d1*d2][d3]
        vector<vector<cfloat>> mat(d1 * d2, vector<cfloat>(d3));
        for (int i = 0; i < d1; i++) {
            for (int j = 0; j < d2; j++) {
                for (int k = 0; k < d3; k++) {
                    mat[i * d2 + j][k] = (*ten3)[i][j][k];
                }
            }
        }
        return mat;
    }

    // Try 4D tensor (rank=4) - flatten to matrix by combining indices pairwise
    if (auto* ten4 = std::get_if<vector<vector<vector<vector<cfloat>>>>>(&result)) {
        if (ten4->empty()) return vector<vector<cfloat>>();
        int d1 = ten4->size();
        int d2 = (*ten4)[0].size();
        int d3 = (*ten4)[0][0].size();
        int d4 = (*ten4)[0][0][0].size();
        // Flatten [d1][d2][d3][d4] to [d1*d2][d3*d4]
        vector<vector<cfloat>> mat(d1 * d2, vector<cfloat>(d3 * d4));
        for (int i = 0; i < d1; i++) {
            for (int j = 0; j < d2; j++) {
                for (int k = 0; k < d3; k++) {
                    for (int l = 0; l < d4; l++) {
                        mat[i * d2 + j][k * d4 + l] = (*ten4)[i][j][k][l];
                    }
                }
            }
        }
        return mat;
    }

    // Return empty matrix on error
    return vector<vector<cfloat>>();
}

vector<float> Field_CM::diag(Vec point, float w) {
    // Get matrix at the specified point
    auto matrix = (*this)(point, w);
    if (matrix.empty()) {
        return vector<float>();
    }

    int N = matrix.size();

    // Copy matrix data into lapack_complex_float array
    vector<lapack_complex_float> A(N * N);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            reinterpret_cast<float(&)[2]>(A[i * N + j])[0] = matrix[i][j].real();
            reinterpret_cast<float(&)[2]>(A[i * N + j])[1] = matrix[i][j].imag();
        }
    }

    // Array to store eigenvalues (real for Hermitian matrices)
    vector<float> eigenvalues(N);

    // Call LAPACK Hermitian eigenvalue solver (only eigenvalues, no eigenvectors)
    int info = LAPACKE_cheev(LAPACK_ROW_MAJOR, 'N', 'U', N, A.data(), N, eigenvalues.data());

    if (info != 0) {
        std::cerr << "Error: LAPACKE_cheev returned " << info << std::endl;
        return vector<float>();
    }

    return eigenvalues;
}

vector<Eigenvector> Field_CM::fulldiag(Vec point, float w) {
    // Get matrix at the specified point
    auto matrix = (*this)(point, w);
    if (matrix.empty()) {
        return vector<Eigenvector>();
    }

    int N = matrix.size();

    // Copy matrix data into lapack_complex_float array
    vector<lapack_complex_float> A(N * N);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            reinterpret_cast<float(&)[2]>(A[i * N + j])[0] = matrix[i][j].real();
            reinterpret_cast<float(&)[2]>(A[i * N + j])[1] = matrix[i][j].imag();
        }
    }

    // Array to store eigenvalues (real for Hermitian matrices)
    vector<float> eigenvalues(N);

    // Call LAPACK Hermitian eigenvalue solver with eigenvectors ('V')
    int info = LAPACKE_cheev(LAPACK_ROW_MAJOR, 'V', 'U', N, A.data(), N, eigenvalues.data());

    if (info != 0) {
        std::cerr << "Error: LAPACKE_cheev returned " << info << std::endl;
        return vector<Eigenvector>();
    }

    // Convert to Eigenvector format
    vector<Eigenvector> eigenvectors(N, Eigenvector(N));
    for (int i = 0; i < N; i++) {
        eigenvectors[i].eigenvalue = eigenvalues[i];
        for (int j = 0; j < N; j++) {
            // Eigenvectors are stored in columns
            eigenvectors[i][j] = reinterpret_cast<float(&)[2]>(A[j * N + i])[0];
        }
    }

    return eigenvectors;
}

vector<vector<vector<cfloat>>> Field_CM::operator()(const vector<Vec>& points, float w) {
    vector<vector<vector<cfloat>>> results(points.size());
    for (size_t i = 0; i < points.size(); i++) {
        results[i] = (*this)(points[i], w);
    }
    return results;
}

vector<vector<cfloat>> Field_CM::operator()(float w) {
    // For matrix fields that are k-independent (or we want k-averaged result),
    // evaluate at the origin point
    Vec origin;
    origin.x = 0;
    origin.y = 0;
    origin.z = 0;
    origin.w = 0;
    origin.area = 0;
    origin.dimension = cmf.data.dimension;
    origin.n = 1;
    auto result = cmf.get_array(origin, w);

    // Try 1D vector (rank=1) - wrap as column matrix
    if (auto* vec = std::get_if<vector<cfloat>>(&result)) {
        vector<vector<cfloat>> mat(vec->size(), vector<cfloat>(1));
        for (size_t i = 0; i < vec->size(); i++) {
            mat[i][0] = (*vec)[i];
        }
        return mat;
    }

    // Try 2D matrix (rank=2)
    if (auto* mat = std::get_if<vector<vector<cfloat>>>(&result)) {
        return *mat;
    }

    // Try 3D tensor (rank=3) - flatten to matrix by combining first two indices
    if (auto* ten3 = std::get_if<vector<vector<vector<cfloat>>>>(&result)) {
        if (ten3->empty()) return vector<vector<cfloat>>();
        int d1 = ten3->size();
        int d2 = (*ten3)[0].size();
        int d3 = (*ten3)[0][0].size();
        vector<vector<cfloat>> mat(d1 * d2, vector<cfloat>(d3));
        for (int i = 0; i < d1; i++) {
            for (int j = 0; j < d2; j++) {
                for (int k = 0; k < d3; k++) {
                    mat[i * d2 + j][k] = (*ten3)[i][j][k];
                }
            }
        }
        return mat;
    }

    // Try 4D tensor (rank=4) - flatten to matrix
    if (auto* ten4 = std::get_if<vector<vector<vector<vector<cfloat>>>>>(&result)) {
        if (ten4->empty()) return vector<vector<cfloat>>();
        int d1 = ten4->size();
        int d2 = (*ten4)[0].size();
        int d3 = (*ten4)[0][0].size();
        int d4 = (*ten4)[0][0][0].size();
        vector<vector<cfloat>> mat(d1 * d2, vector<cfloat>(d3 * d4));
        for (int i = 0; i < d1; i++) {
            for (int j = 0; j < d2; j++) {
                for (int k = 0; k < d3; k++) {
                    for (int l = 0; l < d4; l++) {
                        mat[i * d2 + j][k * d4 + l] = (*ten4)[i][j][k][l];
                    }
                }
            }
        }
        return mat;
    }

    return vector<vector<cfloat>>();
}

vector<vector<vector<cfloat>>> Field_CM::operator()(const vector<float>& w_points) {
    vector<vector<vector<cfloat>>> results(w_points.size());
    for (size_t i = 0; i < w_points.size(); i++) {
        results[i] = (*this)(w_points[i]);
    }
    return results;
}

BaseData* Field_CM::get_data() {
    return &cmf.data;
}

// Field_RM implementation (Real Matrix)
Field_RM::Field_RM()
    : cmf(vector<vector<vector<cfloat>>>(), false, false, true, {}, {}, {}, {3, 3}) {}

Field_RM::Field_RM(const BaseData::DataVariant& data,
                   const vector<int>& inds,
                   const vector<int>& mesh,
                   const vector<vector<float>>& domain,
                   const vector<float>& w_points)
    : cmf(data, false, false, true, mesh, domain, w_points, inds) {}

Field_RM::Field_RM(FieldImpl f) : cmf(f) {}

Field_RM::Field_RM(const string& filename) : cmf(filename) {}

Field_RM& Field_RM::operator=(const Field_RM& other) {
    if (this != &other) {
        cmf = other.cmf;
    }
    return *this;
}

void Field_RM::save(const string& filename) {
    cmf.save(filename);
}

vector<vector<float>> Field_RM::operator()(Vec point, float w) {
    auto result = cmf.get_array(point, w);
    if (auto* mat = std::get_if<vector<vector<float>>>(&result)) {
        return *mat;
    }
    // Return empty matrix on error
    return vector<vector<float>>();
}

vector<float> Field_RM::diag(Vec point, float w) {
    // Get matrix at the specified point
    auto matrix = (*this)(point, w);
    if (matrix.empty()) {
        return vector<float>();
    }

    int N = matrix.size();

    // Copy matrix data into a flat array (LAPACK expects row-major format)
    vector<float> A(N * N);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            A[i * N + j] = matrix[i][j];
        }
    }

    // Array to store eigenvalues
    vector<float> eigenvalues(N);

    // Call LAPACK symmetric eigenvalue solver (only eigenvalues, no eigenvectors)
    int info = LAPACKE_ssyev(LAPACK_ROW_MAJOR, 'N', 'U', N, A.data(), N, eigenvalues.data());

    if (info != 0) {
        std::cerr << "Error: LAPACKE_ssyev returned " << info << std::endl;
        return vector<float>();
    }

    return eigenvalues;
}

vector<Eigenvector> Field_RM::fulldiag(Vec point, float w) {
    // Get matrix at the specified point
    auto matrix = (*this)(point, w);
    if (matrix.empty()) {
        return vector<Eigenvector>();
    }

    int N = matrix.size();

    // Copy matrix data into a flat array (LAPACK expects row-major format)
    vector<float> A(N * N);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            A[i * N + j] = matrix[i][j];
        }
    }

    // Array to store eigenvalues
    vector<float> eigenvalues(N);

    // Call LAPACK symmetric eigenvalue solver with eigenvectors ('V')
    int info = LAPACKE_ssyev(LAPACK_ROW_MAJOR, 'V', 'U', N, A.data(), N, eigenvalues.data());

    if (info != 0) {
        std::cerr << "Error: LAPACKE_ssyev returned " << info << std::endl;
        return vector<Eigenvector>();
    }

    // Convert to Eigenvector format
    // After ssyev, A contains eigenvectors in columns
    vector<Eigenvector> eigenvectors(N, Eigenvector(N));
    for (int i = 0; i < N; i++) {
        eigenvectors[i].eigenvalue = eigenvalues[i];
        for (int j = 0; j < N; j++) {
            // Eigenvectors are stored in columns
            eigenvectors[i][j] = A[j * N + i];
        }
    }

    return eigenvectors;
}

vector<vector<vector<float>>> Field_RM::operator()(const vector<Vec>& points, float w) {
    vector<vector<vector<float>>> results(points.size());
    for (size_t i = 0; i < points.size(); i++) {
        results[i] = (*this)(points[i], w);
    }
    return results;
}

vector<vector<float>> Field_RM::operator()(float w) {
    // For matrix fields that are k-independent (or we want k-averaged result),
    // evaluate at the origin point
    Vec origin;
    origin.x = 0;
    origin.y = 0;
    origin.z = 0;
    origin.w = 0;
    origin.area = 0;
    origin.dimension = cmf.data.dimension;
    origin.n = 1;
    auto result = cmf.get_array(origin, w);

    // Try 1D vector (rank=1) - wrap as column matrix
    if (auto* vec = std::get_if<vector<cfloat>>(&result)) {
        vector<vector<float>> mat(vec->size(), vector<float>(1));
        for (size_t i = 0; i < vec->size(); i++) {
            mat[i][0] = (*vec)[i].real();
        }
        return mat;
    }

    // Try 2D matrix (rank=2)
    if (auto* mat = std::get_if<vector<vector<cfloat>>>(&result)) {
        vector<vector<float>> real_mat(mat->size(), vector<float>((*mat)[0].size()));
        for (size_t i = 0; i < mat->size(); i++) {
            for (size_t j = 0; j < (*mat)[i].size(); j++) {
                real_mat[i][j] = (*mat)[i][j].real();
            }
        }
        return real_mat;
    }

    // Try 3D tensor (rank=3) - flatten to matrix by combining first two indices
    if (auto* ten3 = std::get_if<vector<vector<vector<cfloat>>>>(&result)) {
        if (ten3->empty()) return vector<vector<float>>();
        int d1 = ten3->size();
        int d2 = (*ten3)[0].size();
        int d3 = (*ten3)[0][0].size();
        vector<vector<float>> mat(d1 * d2, vector<float>(d3));
        for (int i = 0; i < d1; i++) {
            for (int j = 0; j < d2; j++) {
                for (int k = 0; k < d3; k++) {
                    mat[i * d2 + j][k] = (*ten3)[i][j][k].real();
                }
            }
        }
        return mat;
    }

    // Try 4D tensor (rank=4) - flatten to matrix
    if (auto* ten4 = std::get_if<vector<vector<vector<vector<cfloat>>>>>(&result)) {
        if (ten4->empty()) return vector<vector<float>>();
        int d1 = ten4->size();
        int d2 = (*ten4)[0].size();
        int d3 = (*ten4)[0][0].size();
        int d4 = (*ten4)[0][0][0].size();
        vector<vector<float>> mat(d1 * d2, vector<float>(d3 * d4));
        for (int i = 0; i < d1; i++) {
            for (int j = 0; j < d2; j++) {
                for (int k = 0; k < d3; k++) {
                    for (int l = 0; l < d4; l++) {
                        mat[i * d2 + j][k * d4 + l] = (*ten4)[i][j][k][l].real();
                    }
                }
            }
        }
        return mat;
    }

    return vector<vector<float>>();
}

vector<vector<vector<float>>> Field_RM::operator()(const vector<float>& w_points) {
    vector<vector<vector<float>>> results(w_points.size());
    for (size_t i = 0; i < w_points.size(); i++) {
        results[i] = (*this)(w_points[i]);
    }
    return results;
}

BaseData* Field_RM::get_data() {
    return &cmf.data;
}

// UnifiedField implementation
Field::Field() {
    is_complex = false;
    is_vector = false;
    is_matrix = false;
    field_c = nullptr;
    field_r = nullptr;
    field_cm = nullptr;
    field_rm = nullptr;
    default_plot_type = "";
    title = "";
    x_label = "";
    y_label = "";
}

void Field::generate_plot_labels(const string& filename) {
    // Reset labels
    title = "";
    x_label = "";
    y_label = "";
    default_plot_type = "line";

    if (filename.empty()) return;

    // Extract directory and basename
    size_t last_slash = filename.find_last_of("/\\");
    string full_dir = (last_slash != string::npos) ? filename.substr(0, last_slash) : "";
    string base_name = (last_slash != string::npos) ? filename.substr(last_slash + 1) : filename;

    // Extract only the immediate directory name (not full path)
    string dir_part = "";
    if (!full_dir.empty()) {
        size_t second_last_slash = full_dir.find_last_of("/\\");
        dir_part = (second_last_slash != string::npos) ? full_dir.substr(second_last_slash + 1) : full_dir;
    }

    // Remove .h5 extension
    size_t dot_pos = base_name.find_last_of('.');
    if (dot_pos != string::npos) {
        base_name = base_name.substr(0, dot_pos);
    }

    // Parse pattern: prefix_field_vars
    size_t last_underscore = base_name.find_last_of('_');
    if (last_underscore == string::npos) return;

    string vars = base_name.substr(last_underscore + 1);
    string prefix_and_field = base_name.substr(0, last_underscore);

    // Extract field name (after last underscore before vars)
    size_t second_last_underscore = prefix_and_field.find_last_of('_');
    string field_name;
    string prefix;

    if (second_last_underscore != string::npos) {
        field_name = prefix_and_field.substr(second_last_underscore + 1);
        prefix = prefix_and_field.substr(0, second_last_underscore);
    } else {
        field_name = prefix_and_field;
        prefix = "";
    }

    // Map field names to LaTeX
    string field_latex = "";
    if (field_name == "sigma" || field_name == "Sigma") {
        field_latex = "\\Sigma";
    } else if (field_name == "chi" || field_name == "Chi") {
        field_latex = "\\chi";
    } else if (field_name == "vertex") {
        field_latex = "V";
    } else if (field_name == "G") {
        field_latex = "G";
    } else if (field_name == "Delta" || field_name == "delta") {
        field_latex = "\\Delta";
    } else if (field_name == "hamiltonian" || field_name == "H") {
        field_latex = "H";
    } else {
        field_latex = field_name;
    }

    // Parse variables and create argument list
    string arg_list = "";
    bool has_w = false, has_k = false, has_q = false;

    for (char c : vars) {
        if (c == 'w') {
            has_w = true;
            if (!arg_list.empty()) arg_list += ",";
            arg_list += "\\omega";
        } else if (c == 'k') {
            has_k = true;
            if (!arg_list.empty()) arg_list += ",";
            arg_list += "k";
        } else if (c == 'q') {
            has_q = true;
            if (!arg_list.empty()) arg_list += ",";
            arg_list += "q";
        } else if (c == 't') {
            if (!arg_list.empty()) arg_list += ",";
            arg_list += "\\tau";
        } else if (c == 'r') {
            if (!arg_list.empty()) arg_list += ",";
            arg_list += "r";
        }
    }

    // Create y_label
    if (!field_latex.empty() && !arg_list.empty()) {
        y_label = "$" + field_latex + "(" + arg_list + ")$";
    } else if (!field_latex.empty()) {
        y_label = "$" + field_latex + "$";
    }

    // Create x_label (prioritize: w > k > q)
    if (has_w) {
        x_label = "$\\omega$";
    } else if (has_k) {
        x_label = "$k$";
    } else if (has_q) {
        x_label = "$q$";
    }

    // Create title
    if (!y_label.empty()) {
        if (!dir_part.empty() && !prefix.empty()) {
            title = y_label + " for " + dir_part + "/" + prefix;
        } else if (!dir_part.empty()) {
            title = y_label + " for " + dir_part;
        } else if (!prefix.empty()) {
            title = y_label + " for " + prefix;
        } else {
            title = y_label;
        }
    }

    // Set default plot type based on variables
    if (has_w && has_k) {
        default_plot_type = "heatmap";
    } else if (has_k || has_w) {
        default_plot_type = "line";
    } else {
        default_plot_type = "scatter";
    }
}

Field::Field(const string& filename) {
    // Load metadata from file to determine type
    BaseData base = load_data_from_hdf5(filename);
    is_complex = base.is_complex;
    is_vector = base.is_vector;
    is_matrix = base.is_matrix;

    field_c = nullptr;
    field_r = nullptr;
    field_cm = nullptr;
    field_rm = nullptr;

    // Create appropriate field type
    if (is_matrix && is_complex) {
        field_cm = new Field_CM(filename);
    } else if (is_matrix && !is_complex) {
        field_rm = new Field_RM(filename);
    } else if (!is_matrix && is_complex) {
        field_c = new Field_C(filename);
    } else {
        field_r = new Field_R(filename);
    }

    // Generate plot labels from filename
    generate_plot_labels(filename);
}

Field::~Field() {
    if (field_c) delete field_c;
    if (field_r) delete field_r;
    if (field_cm) delete field_cm;
    if (field_rm) delete field_rm;
}

void Field::save(const string& filename) {
    if (field_c) field_c->save(filename);
    else if (field_r) field_r->save(filename);
    else if (field_cm) field_cm->save(filename);
    else if (field_rm) field_rm->save(filename);
}

// Scalar complex operators
complex<float> Field::operator_scalar_complex(Vec point, float w) {
    if (!field_c) throw runtime_error("Field is not complex scalar");
    return (*field_c)(point, w);
}

complex<float> Field::operator_scalar_complex(float w) {
    if (!field_c) throw runtime_error("Field is not complex scalar");
    return (*field_c)(w);
}

vector<complex<float>> Field::operator_scalar_complex(const vector<Vec>& points, float w) {
    if (!field_c) throw runtime_error("Field is not complex scalar");
    return (*field_c)(points, w);
}

vector<complex<float>> Field::operator_scalar_complex(const vector<float>& w_points) {
    if (!field_c) throw runtime_error("Field is not complex scalar");
    return (*field_c)(w_points);
}

// Scalar real operators
float Field::operator_scalar_real(Vec point, float w) {
    if (!field_r) throw runtime_error("Field is not real scalar");
    return (*field_r)(point, w);
}

float Field::operator_scalar_real(float w) {
    if (!field_r) throw runtime_error("Field is not real scalar");
    return (*field_r)(w);
}

vector<float> Field::operator_scalar_real(const vector<Vec>& points, float w) {
    if (!field_r) throw runtime_error("Field is not real scalar");
    return (*field_r)(points, w);
}

vector<float> Field::operator_scalar_real(const vector<float>& w_points) {
    if (!field_r) throw runtime_error("Field is not real scalar");
    return (*field_r)(w_points);
}

// Matrix complex operators
vector<vector<complex<float>>> Field::operator_matrix_complex(Vec point, float w) {
    if (!field_cm) throw runtime_error("Field is not complex matrix");
    return (*field_cm)(point, w);
}

vector<vector<vector<complex<float>>>> Field::operator_matrix_complex(const vector<Vec>& points, float w) {
    if (!field_cm) throw runtime_error("Field is not complex matrix");
    return (*field_cm)(points, w);
}

// Matrix real operators
vector<vector<float>> Field::operator_matrix_real(Vec point, float w) {
    if (!field_rm) throw runtime_error("Field is not real matrix");
    return (*field_rm)(point, w);
}

vector<vector<vector<float>>> Field::operator_matrix_real(const vector<Vec>& points, float w) {
    if (!field_rm) throw runtime_error("Field is not real matrix");
    return (*field_rm)(points, w);
}

BaseData* Field::get_data() {
    if (field_c) return field_c->get_data();
    if (field_r) return field_r->get_data();
    if (field_cm) return field_cm->get_data();
    if (field_rm) return field_rm->get_data();
    return nullptr;
}
