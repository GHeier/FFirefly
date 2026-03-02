#include "src/objects/CMField/fields.hpp"
#include "src/objects/CMField/base_data.hpp"
#include <filesystem>
#include <iostream>

#include "src/config/load/c_config.h"

using namespace std;

static int mpts = 3;

static cfloat func_quad(Vec p, int dim) {
    float val = 0.0;
    for (int i = 0; i < dim; i++)
        val += cos(p(i)) * cos(p(i));
    return cfloat(val, 0);
}

static float get_pnt(int i, int pnts) {
    return -M_PI + 2 * M_PI * i / pnts;
}

static float get_pnt_fft(int i, int pnts) {
    return 2 * M_PI * i / pnts;
}

static Vec get_vec(int i, int j, int k, int pnts, bool fft = false) {
    if (fft)
    return Vec(
            get_pnt_fft(i, pnts),
            get_pnt_fft(j, pnts),
            get_pnt_fft(k, pnts)
            );
    return Vec(
            get_pnt(i, pnts),
            get_pnt(j, pnts),
            get_pnt(k, pnts)
            );
}

vector<cfloat> create_data(int dim, int pnts, bool fft = false) {
    vector<cfloat> values;
    for (int i = 0; i < pnts; i++) {
        for (int j = 0; j < pnts; j++) {
            Vec point = get_vec(i, j, 0, pnts, fft);
            point.dimension = 2;
            cfloat base = func_quad(point, dim);
            values.push_back(base);
            //cout << point << " : " << base.real() << endl;
        }
    }
    return values;
}

static bool field_r_2d_k() {
    vector<int> mesh = {mpts, mpts};
    vector<vector<float>> domain = {{2 * M_PI, 0}, {0, 2 * M_PI}};
    vector<cfloat> data = create_data(2, mpts);

    //BaseData::DataVariant sd = data;

    Field_C field(data, mesh, domain);

    Vec v(1.0472, 1.0472);
    float result = field(v).real();
    float expected = 0.5;

    save_data("testfield.h5", data, {}, mesh, domain);
    Field_C sfield("testfield.h5");
    result = sfield(v).real();

    field.save("testfield.h5");
    Field_C fsfield("testfield.h5");
    result = fsfield(v).real();

    return fabs(result - expected) < 1e-5;
}

static bool field_r_2d_k_fft() {
    vector<int> mesh = {mpts, mpts};
    vector<vector<float>> domain = {{2 * M_PI, 0}, {0, 2 * M_PI}};
    vector<cfloat> data = create_data(2, mpts, true);

    //BaseData::DataVariant sd = data;

    Field_C field(data, mesh, domain, {}, true);

    Vec v(1.0472, 1.0472);
    float result = field(v).real();
    float expected = 0.5;

    save_data("testfield.h5", data, {}, mesh, domain);
    Field_C sfield("testfield.h5");
    result = sfield(v).real();

    field.save("testfield.h5");
    Field_C fsfield("testfield.h5");
    result = fsfield(v).real();

    return fabs(result - expected) < 1e-5;
}


bool quad_field_tests() {
    int num_tests = 2;
    bool all_tests[num_tests] = {
        field_r_2d_k(),
        field_r_2d_k_fft(),
    };
    filesystem::remove("testfield.h5");
    return print_test_results(all_tests, num_tests, "Quadratic Field tests");
}
