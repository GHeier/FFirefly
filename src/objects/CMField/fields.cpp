#include "fields.hpp"
#include "field.hpp"
#include "../vec.hpp"

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
