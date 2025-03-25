#include "../cmdata.hpp"
#include "../../vec.hpp"
#include "all.hpp"

#include "../../../config/load/c_config.h"

using namespace std;

bool CMData_tests() {
    printf("Running CMData tests\n");
    int num_tests = 3;
    bool all_tests[num_tests] = {
        test_1d(),
        test_2d(),
        test_3d(),
    };
    return print_test_results(all_tests, num_tests, "CMData tests");
}

int pts = 3;

complex<Vec> func_d(Vec x) {
    Vec val(x.norm() * x.norm());
    Vec empty;
    complex<Vec> r_val = complex<Vec>(val, empty);
    return r_val;
}

float indf(int i) {
    return 2.0 * i / pts - 1.0 * i / pts;
}

bool test_1d() {
    vector<Vec> points;
    vector<complex<Vec>> values;
    for (int i = 0; i < pts; i++) {
        Vec x(indf(i));
        points.push_back(x);
        values.push_back(func_d(x));
    }
    float test_pnt = points[2](0);
    float test_val = real(values[2])(0);

    CMData data(points, values, 1, false, false, false, false);

    float data_pnt = data.points[2](0);
    float data_val = real(data.values[2])(0);

    if (abs(test_pnt - data_pnt) > 1e-6) return false;
    if (abs(test_val - data_val) > 1e-6) return false;

    save(data, "testfile");
    CMData loaddata("testfile");

    float load_pnt = loaddata.points[2](0);
    float load_val = real(loaddata.values[2])(0);

    if (abs(test_pnt - load_pnt) > 1e-6) return false;
    if (abs(test_val - load_val) > 1e-6) return false;

    return true;
}

bool test_2d() {
    vector<Vec> points;
    vector<complex<Vec>> values;
    for (int i = 0; i < pts; i++) {
        for (int j = 0; j < pts; j++) {
            Vec x(indf(i));
            points.push_back(x);
            values.push_back(func_d(x));
        }
    }
    float test_pnt = points[2](0);
    float test_val = real(values[2])(0);

    CMData data(points, values, 1, false, false, false, false);

    float data_pnt = data.points[2](0);
    float data_val = real(data.values[2])(0);

    if (abs(test_pnt - data_pnt) > 1e-6) return false;
    if (abs(test_val - data_val) > 1e-6) return false;

    save(data, "testfile");
    CMData loaddata("testfile");

    float load_pnt = loaddata.points[2](0);
    float load_val = real(loaddata.values[2])(0);

    if (abs(test_pnt - load_pnt) > 1e-6) return false;
    if (abs(test_val - load_val) > 1e-6) return false;

    return true;
}

bool test_3d() {
    vector<Vec> points;
    vector<complex<Vec>> values;
    for (int i = 0; i < pts; i++) {
        for (int j = 0; j < pts; j++) {
            for (int k = 0; k < pts; k++) {
                Vec x(indf(i));
                points.push_back(x);
                values.push_back(func_d(x));
            }
        }
    }
    float test_pnt = points[2](0);
    float test_val = real(values[2])(0);

    CMData data(points, values, 1, false, false, false, false);

    float data_pnt = data.points[2](0);
    float data_val = real(data.values[2])(0);

    if (abs(test_pnt - data_pnt) > 1e-6) return false;
    if (abs(test_val - data_val) > 1e-6) return false;

    save(data, "testfile");
    CMData loaddata("testfile");

    float load_pnt = loaddata.points[2](0);
    float load_val = real(loaddata.values[2])(0);

    if (abs(test_pnt - load_pnt) > 1e-6) return false;
    if (abs(test_val - load_val) > 1e-6) return false;

    return true;
}
