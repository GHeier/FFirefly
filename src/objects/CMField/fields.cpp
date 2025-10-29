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

Field_C::Field_C(Field f) : cmf(f) {}

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

Field_R::Field_R(Field f) : cmf(f) {}

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
    : cmf(vector<vector<vector<cfloat>>>(), true, false, true, {}, {}, {}, 2, 1) {}

Field_CM::Field_CM(const BaseData::DataVariant& data,
                   int dim_indices,
                   const vector<int>& mesh,
                   const vector<vector<float>>& domain,
                   const vector<float>& w_points)
    : cmf(data, true, false, true, mesh, domain, w_points, 2, dim_indices) {}

Field_CM::Field_CM(Field f) : cmf(f) {}

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
    if (auto* mat = std::get_if<vector<vector<cfloat>>>(&result)) {
        return *mat;
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

BaseData* Field_CM::get_data() {
    return &cmf.data;
}

// Field_RM implementation (Real Matrix)
Field_RM::Field_RM()
    : cmf(vector<vector<vector<cfloat>>>(), false, false, true, {}, {}, {}, 2, 1) {}

Field_RM::Field_RM(const BaseData::DataVariant& data,
                   int dim_indices,
                   const vector<int>& mesh,
                   const vector<vector<float>>& domain,
                   const vector<float>& w_points)
    : cmf(data, false, false, true, mesh, domain, w_points, 2, dim_indices) {}

Field_RM::Field_RM(Field f) : cmf(f) {}

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

BaseData* Field_RM::get_data() {
    return &cmf.data;
}
