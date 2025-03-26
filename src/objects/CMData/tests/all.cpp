#include "../cmdata.hpp"
#include "../../vec.hpp"
#include "all.hpp"

#include "../../../config/load/c_config.h"

using namespace std;

bool CMData_tests() {
    printf("Running CMData tests\n");
    int num_tests = 27;
    bool all_tests[num_tests] = {
        readwrite_1d(),
        readwrite_2d(),
        readwrite_3d(),
        readwrite_1d_complex(),
        readwrite_2d_complex(),
        readwrite_3d_complex(),
        readwrite_0d_with_w(),
        readwrite_1d_with_w(),
        readwrite_2d_with_w(),
        readwrite_3d_with_w(),
        readwrite_1d_vector(),
        readwrite_2d_vector(),
        readwrite_3d_vector(),
        readwrite_1d_complexvector(),
        readwrite_2d_complexvector(),
        readwrite_3d_complexvector(),
        readwrite_0d_complexvector_with_w(),
        readwrite_1d_complexvector_with_w(),
        readwrite_2d_complexvector_with_w(),
        readwrite_3d_complexvector_with_w(),
        readwrite_1d_with_n(),
        readwrite_2d_with_n(),
        readwrite_3d_with_n(),
        readwrite_0d_with_w_with_n(),
        readwrite_1d_with_w_with_n(),
        readwrite_2d_with_w_with_n(),
        readwrite_3d_with_w_with_n(),
    };
    remove("testfile");
    return print_test_results(all_tests, num_tests, "CMData tests");
}

int pts = 3;
int npts = 3;

complex<Vec> func_d(Vec x, int dim) {
    Vec val(x.norm() * x.norm(), x.norm() * x.norm(), x.norm() * x.norm());
    val.dimension = dim;
    complex<Vec> r_val = complex<Vec>(val, val);
    return r_val;
}

float indf(int i) {
    return 2.0 * i / pts - 1.0;
}

bool readwrite_1d() {
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

bool readwrite_1d_complex() {
    vector<Vec> points;
    vector<complex<Vec>> values;
    for (int i = 0; i < pts; i++) {
        Vec x(indf(i));
        points.push_back(x);
        values.push_back(func_d(x));
    }
    float test_pnt = points[2](0);
    complex<float> test_val = complex<float>(real(values[2])(0), imag(values[2])(0));

    CMData data(points, values, 1, false, false, true, false);

    float data_pnt = data.points[2](0);
    complex<float> data_val = complex<float>(real(data.values[2])(0), imag(data.values[2])(0));

    if (abs(test_pnt - data_pnt) > 1e-6) return false;
    if (abs(test_val - data_val) > 1e-6) return false;

    save(data, "testfile");
    CMData loaddata("testfile");

    float load_pnt = loaddata.points[2](0);
    complex<float> load_val = complex<float>(real(loaddata.values[2])(0), imag(loaddata.values[2])(0));

    if (abs(test_pnt - load_pnt) > 1e-6) return false;
    if (abs(test_val - load_val) > 1e-6) return false;

    return true;
}

bool readwrite_0d_with_w() {
    vector<Vec> points;
    vector<complex<Vec>> values;
    for (int i = 0; i < pts; i++) {
        Vec x(indf(i));
        points.push_back(x);
        values.push_back(func_d(x));
    }
    float test_pnt = points[2](0);
    float test_val = real(values[2])(0);

    CMData data(points, values, 0, true, false, false, false);

    float data_pnt = data.points[2](0);
    float data_val = real(data.values[2])(0);
    float data_wpnt = data.w_points[2];

    if (abs(test_pnt - data_pnt) > 1e-6) return false;
    if (abs(test_val - data_val) > 1e-6) return false;
    if (abs(test_pnt - data_wpnt) > 1e-6) return false;

    save(data, "testfile");
    CMData loaddata("testfile");

    float load_pnt = loaddata.points[2](0);
    float load_val = real(loaddata.values[2])(0);
    float load_wpnt = loaddata.w_points[2];

    if (abs(test_pnt - load_pnt) > 1e-6) return false;
    if (abs(test_val - load_val) > 1e-6) return false;
    if (abs(test_pnt - load_wpnt) > 1e-6) return false;

    return true;
}

bool readwrite_1d_with_w() {
    vector<Vec> points;
    vector<complex<Vec>> values;
    for (int i = 0; i < pts; i++) {
        for (int j = 0; j < pts; j++) {
            Vec x(indf(i), indf(j));
            points.push_back(x);
            values.push_back(func_d(x));
        }
    }
    float test_pnt = points[2](0);
    float test_val = real(values[2])(0);
    float test_wpnt = points[2](1);

    CMData data(points, values, 1, true, false, false, false);

    float data_pnt = data.points[2](0);
    float data_val = real(data.values[2])(0);
    float data_wpnt = data.w_points[2];

    if (abs(test_pnt - data_pnt) > 1e-6) return false;
    if (abs(test_val - data_val) > 1e-6) return false;
    if (abs(test_wpnt - data_wpnt) > 1e-6) return false;

    save(data, "testfile");
    CMData loaddata("testfile");

    float load_pnt = loaddata.points[2](0);
    float load_val = real(loaddata.values[2])(0);
    float load_wpnt = loaddata.w_points[2];

    if (abs(test_pnt - load_pnt) > 1e-6) return false;
    if (abs(test_val - load_val) > 1e-6) return false;
    if (abs(test_wpnt - load_wpnt) > 1e-6) return false;

    return true;
}

bool readwrite_2d_with_w() {
    vector<Vec> points;
    vector<complex<Vec>> values;
    for (int i = 0; i < pts; i++) {
        for (int j = 0; j < pts; j++) {
            for (int k = 0; k < pts; k++) {
                Vec x(indf(i), indf(j), indf(k));
                points.push_back(x);
                values.push_back(func_d(x));
            }
        }
    }
    float test_pnt = points[2](0);
    float test_val = real(values[2])(0);
    float test_wpnt = points[2](2);

    CMData data(points, values, 2, true, false, false, false);

    float data_pnt = data.points[2](0);
    float data_val = real(data.values[2])(0);
    float data_wpnt = data.w_points[2];

    if (abs(test_pnt - data_pnt) > 1e-6) return false;
    if (abs(test_val - data_val) > 1e-6) return false;
    if (abs(test_wpnt - data_wpnt) > 1e-6) return false;

    save(data, "testfile");
    CMData loaddata("testfile");

    float load_pnt = loaddata.points[2](0);
    float load_val = real(loaddata.values[2])(0);
    float load_wpnt = loaddata.w_points[2];

    if (abs(test_pnt - load_pnt) > 1e-6) return false;
    if (abs(test_val - load_val) > 1e-6) return false;
    if (abs(test_wpnt - load_wpnt) > 1e-6) return false;

    return true;
}

bool readwrite_3d_with_w() {
    vector<Vec> points;
    vector<complex<Vec>> values;
    for (int i = 0; i < pts; i++) {
        for (int j = 0; j < pts; j++) {
            for (int k = 0; k < pts; k++) {
                for (int l = 0; l < pts; l++) {
                    Vec x(indf(i), indf(j), indf(k), indf(l));
                    points.push_back(x);
                    values.push_back(func_d(x));
                }
            }
        }
    }
    float test_pnt = points[2](0);
    float test_val = real(values[2])(0);
    float test_wpnt = points[2](3);

    CMData data(points, values, 3, true, false, false, false);

    float data_pnt = data.points[2](0);
    float data_val = real(data.values[2])(0);
    float data_wpnt = data.w_points[2];

    if (abs(test_pnt - data_pnt) > 1e-6) return false;
    if (abs(test_val - data_val) > 1e-6) return false;
    if (abs(test_wpnt - data_wpnt) > 1e-6) return false;

    save(data, "testfile");
    CMData loaddata("testfile");

    float load_pnt = loaddata.points[2](0);
    float load_val = real(loaddata.values[2])(0);
    float load_wpnt = loaddata.w_points[2];

    if (abs(test_pnt - load_pnt) > 1e-6) return false;
    if (abs(test_val - load_val) > 1e-6) return false;
    if (abs(test_wpnt - load_wpnt) > 1e-6) return false;

    return true;
}

bool readwrite_2d() {
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

    CMData data(points, values, 2, false, false, false, false);

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

bool readwrite_2d_complex() {
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

    CMData data(points, values, 2, false, false, true, false);

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

bool readwrite_3d() {
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

    CMData data(points, values, 3, false, false, false, false);

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

bool readwrite_3d_complex() {
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

    CMData data(points, values, 3, false, false, true, false);

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

bool readwrite_1d_vector() {
    vector<Vec> points;
    vector<complex<Vec>> values;
    for (int i = 0; i < pts; i++) {
        Vec x(indf(i));
        x.dimension = 1;
        points.push_back(x);
        values.push_back(func_d(x, x.dimension));
    }
    float test_pnt = points[2](0);
    Vec test_val = real(values[2]);

    CMData data(points, values, 1, false, false, false, true);

    float data_pnt = data.points[2](0);
    Vec data_val = real(data.values[2]);

    if (abs(test_pnt - data_pnt) > 1e-6) return false;
    if ((test_val - data_val).norm() > 1e-6) return false;

    save(data, "testfile");
    CMData loaddata("testfile");

    float load_pnt = loaddata.points[2](0);
    Vec load_val = real(loaddata.values[2]);

    if (abs(test_pnt - load_pnt) > 1e-6) return false;
    if ((test_val - load_val).norm() > 1e-6) return false;

    return true;
}

bool readwrite_2d_vector() {
    vector<Vec> points;
    vector<complex<Vec>> values;
    for (int i = 0; i < pts; i++) {
        for (int j = 0; j < pts; j++) {
            Vec x(indf(i), indf(j));
            x.dimension = 2;
            points.push_back(x);
            values.push_back(func_d(x, x.dimension));
        }
    }
    float test_pnt = points[2](0);
    Vec test_val = real(values[2]);

    CMData data(points, values, 2, false, false, false, true);

    float data_pnt = data.points[2](0);
    Vec data_val = real(data.values[2]);

    if (abs(test_pnt - data_pnt) > 1e-6) return false;
    if ((test_val - data_val).norm() > 1e-6) return false;

    save(data, "testfile");
    CMData loaddata("testfile");

    float load_pnt = loaddata.points[2](0);
    Vec load_val = real(loaddata.values[2]);

    if (abs(test_pnt - load_pnt) > 1e-6) return false;
    if ((test_val - load_val).norm() > 1e-6) return false;

    return true;
}

bool readwrite_3d_vector() {
    vector<Vec> points;
    vector<complex<Vec>> values;
    for (int i = 0; i < pts; i++) {
        for (int j = 0; j < pts; j++) {
            for (int k = 0; k < pts; k++) {
                Vec x(indf(i), indf(j), indf(k));
                x.dimension = 3;
                points.push_back(x);
                values.push_back(func_d(x, x.dimension));
            }
        }
    }
    float test_pnt = points[2](0);
    Vec test_val = real(values[2]);

    CMData data(points, values, 3, false, false, false, true);

    float data_pnt = data.points[2](0);
    Vec data_val = real(data.values[2]);

    if (abs(test_pnt - data_pnt) > 1e-6) return false;
    if ((test_val - data_val).norm() > 1e-6) return false;

    save(data, "testfile");
    CMData loaddata("testfile");

    float load_pnt = loaddata.points[2](0);
    Vec load_val = real(loaddata.values[2]);

    if (abs(test_pnt - load_pnt) > 1e-6) return false;
    if ((test_val - load_val).norm() > 1e-6) return false;

    return true;
}

bool readwrite_1d_complexvector() {
    vector<Vec> points;
    vector<complex<Vec>> values;
    for (int i = 0; i < pts; i++) {
        Vec x(indf(i));
        points.push_back(x);
        values.push_back(func_d(x, 1));
    }
    float test_pnt = points[2](0);
    complex<Vec> test_val = values[2];

    CMData data(points, values, 1, false, false, true, true);

    float data_pnt = data.points[2](0);
    complex<Vec> data_val = data.values[2];

    if (abs(test_pnt - data_pnt) > 1e-6) return false;
    if ((real(test_val) - real(data_val)).norm() > 1e-6 and (imag(test_val) - imag(data_val)).norm() > 1e-6) return false;

    save(data, "testfile");
    CMData loaddata("testfile");

    float load_pnt = loaddata.points[2](0);
    complex<Vec> load_val = loaddata.values[2];

    if (abs(test_pnt - load_pnt) > 1e-6) return false;
    if ((real(test_val) - real(load_val)).norm() > 1e-6 and (imag(test_val) - imag(load_val)).norm() > 1e-6) return false;

    return true;
}

bool readwrite_2d_complexvector() {
    vector<Vec> points;
    vector<complex<Vec>> values;
    for (int i = 0; i < pts; i++) {
        for (int j = 0; j < pts; j++) {
            Vec x(indf(i), indf(j));
            points.push_back(x);
            values.push_back(func_d(x, 2));
        }
    }
    float test_pnt = points[2](0);
    complex<Vec> test_val = values[2];

    CMData data(points, values, 2, false, false, true, true);

    float data_pnt = data.points[2](0);
    complex<Vec> data_val = data.values[2];

    if (abs(test_pnt - data_pnt) > 1e-6) return false;
    if ((real(test_val) - real(data_val)).norm() > 1e-6 and (imag(test_val) - imag(data_val)).norm() > 1e-6) return false;

    save(data, "testfile");
    CMData loaddata("testfile");

    float load_pnt = loaddata.points[2](0);
    complex<Vec> load_val = loaddata.values[2];

    if (abs(test_pnt - load_pnt) > 1e-6) return false;
    if ((real(test_val) - real(load_val)).norm() > 1e-6 and (imag(test_val) - imag(load_val)).norm() > 1e-6) return false;

    return true;
}

bool readwrite_3d_complexvector() {
    vector<Vec> points;
    vector<complex<Vec>> values;
    for (int i = 0; i < pts; i++) {
        for (int j = 0; j < pts; j++) {
            for (int k = 0; k < pts; k++) {
                Vec x(indf(i), indf(j), indf(k));
                points.push_back(x);
                values.push_back(func_d(x, 3));
            }
        }
    }
    float test_pnt = points[2](0);
    complex<Vec> test_val = values[2];

    CMData data(points, values, 3, false, false, true, true);

    float data_pnt = data.points[2](0);
    complex<Vec> data_val = data.values[2];

    if (abs(test_pnt - data_pnt) > 1e-6) return false;
    if ((real(test_val) - real(data_val)).norm() > 1e-6 and (imag(test_val) - imag(data_val)).norm() > 1e-6) return false;

    save(data, "testfile");
    CMData loaddata("testfile");

    float load_pnt = loaddata.points[2](0);
    complex<Vec> load_val = loaddata.values[2];

    if (abs(test_pnt - load_pnt) > 1e-6) return false;
    if ((real(test_val) - real(load_val)).norm() > 1e-6 and (imag(test_val) - imag(load_val)).norm() > 1e-6) return false;

    return true;
}

bool readwrite_0d_complexvector_with_w() {
    vector<Vec> points;
    vector<complex<Vec>> values;
    for (int i = 0; i < pts; i++) {
        Vec x(indf(i));
        points.push_back(x);
        values.push_back(func_d(x, 1));
    }
    float test_pnt = points[2](0);
    complex<Vec> test_val = values[2];

    CMData data(points, values, 0, true, false, true, true);

    float data_pnt = data.points[2](0);
    complex<Vec> data_val = data.values[2];

    if (abs(test_pnt - data_pnt) > 1e-6) return false;
    if ((real(test_val) - real(data_val)).norm() > 1e-6 and (imag(test_val) - imag(data_val)).norm() > 1e-6) return false;

    save(data, "testfile");
    CMData loaddata("testfile");

    float load_pnt = loaddata.points[2](0);
    complex<Vec> load_val = loaddata.values[2];

    if (abs(test_pnt - load_pnt) > 1e-6) return false;
    if ((real(test_val) - real(load_val)).norm() > 1e-6 and (imag(test_val) - imag(load_val)).norm() > 1e-6) return false;

    return true;
}

bool readwrite_1d_complexvector_with_w() {
    vector<Vec> points;
    vector<complex<Vec>> values;
    for (int i = 0; i < pts; i++) {
        Vec x(indf(i));
        points.push_back(x);
        values.push_back(func_d(x, 1));
    }
    float test_pnt = points[2](0);
    complex<Vec> test_val = values[2];

    CMData data(points, values, 1, true, false, true, true);

    float data_pnt = data.points[2](0);
    complex<Vec> data_val = data.values[2];

    if (abs(test_pnt - data_pnt) > 1e-6) return false;
    if ((real(test_val) - real(data_val)).norm() > 1e-6 and (imag(test_val) - imag(data_val)).norm() > 1e-6) return false;

    save(data, "testfile");
    CMData loaddata("testfile");

    float load_pnt = loaddata.points[2](0);
    complex<Vec> load_val = loaddata.values[2];

    if (abs(test_pnt - load_pnt) > 1e-6) return false;
    if ((real(test_val) - real(load_val)).norm() > 1e-6 and (imag(test_val) - imag(load_val)).norm() > 1e-6) return false;

    return true;
}

bool readwrite_2d_complexvector_with_w() {
    vector<Vec> points;
    vector<complex<Vec>> values;
    for (int i = 0; i < pts; i++) {
        for (int j = 0; j < pts; j++) {
            Vec x(indf(i), indf(j));
            points.push_back(x);
            values.push_back(func_d(x, 2));
        }
    }
    float test_pnt = points[2](0);
    complex<Vec> test_val = values[2];

    CMData data(points, values, 2, true, false, true, true);

    float data_pnt = data.points[2](0);
    complex<Vec> data_val = data.values[2];

    if (abs(test_pnt - data_pnt) > 1e-6) return false;
    if ((real(test_val) - real(data_val)).norm() > 1e-6 and (imag(test_val) - imag(data_val)).norm() > 1e-6) return false;

    save(data, "testfile");
    CMData loaddata("testfile");

    float load_pnt = loaddata.points[2](0);
    complex<Vec> load_val = loaddata.values[2];

    if (abs(test_pnt - load_pnt) > 1e-6) return false;
    if ((real(test_val) - real(load_val)).norm() > 1e-6 and (imag(test_val) - imag(load_val)).norm() > 1e-6) return false;

    return true;
}

bool readwrite_3d_complexvector_with_w() {
    vector<Vec> points;
    vector<complex<Vec>> values;
    for (int i = 0; i < pts; i++) {
        for (int j = 0; j < pts; j++) {
            for (int k = 0; k < pts; k++) {
                Vec x(indf(i), indf(j), indf(k));
                points.push_back(x);
                values.push_back(func_d(x, 3));
            }
        }
    }
    float test_pnt = points[2](0);
    complex<Vec> test_val = values[2];

    CMData data(points, values, 3, true, false, true, true);

    float data_pnt = data.points[2](0);
    complex<Vec> data_val = data.values[2];

    if (abs(test_pnt - data_pnt) > 1e-6) return false;
    if ((real(test_val) - real(data_val)).norm() > 1e-6 and (imag(test_val) - imag(data_val)).norm() > 1e-6) return false;

    save(data, "testfile");
    CMData loaddata("testfile");

    float load_pnt = loaddata.points[2](0);
    complex<Vec> load_val = loaddata.values[2];

    if (abs(test_pnt - load_pnt) > 1e-6) return false;
    if ((real(test_val) - real(load_val)).norm() > 1e-6 and (imag(test_val) - imag(load_val)).norm() > 1e-6) return false;

    return true;
}
bool readwrite_1d_with_n() {
    vector<Vec> points;
    vector<complex<Vec>> values;
    vector<int> ninds;
    for (int n = 1; n <= npts; n++) {
        for (int i = 0; i < pts; i++) {
            Vec x(indf(i));
            x.dimension = 1;
            x.n = n;
            points.push_back(x);
            values.push_back(func_d(x));
        }
        ninds.push_back((n-1) * pts);
    }
    Vec test_pnt = points[2];
    float test_val = real(values[2])(0);

    CMData data(points, values, 1, false, true, false, false);

    if (ninds != data.n_inds) return false;

    Vec data_pnt = data.points[2];
    float data_val = real(data.values[2])(0);

    if ((test_pnt - data_pnt).norm() > 1e-6) return false;
    if (abs(test_val - data_val) > 1e-6) return false;

    save(data, "testfile");
    CMData loaddata("testfile");

    Vec load_pnt = loaddata.points[2];
    float load_val = real(loaddata.values[2])(0);

    if ((test_pnt - load_pnt).norm() > 1e-6) return false;
    if (abs(test_val - load_val) > 1e-6) return false;

    return true;
}

bool readwrite_2d_with_n() {
    vector<Vec> points;
    vector<complex<Vec>> values;
    vector<int> ninds;
    for (int n = 1; n <= npts; n++) {
        for (int i = 0; i < pts; i++) {
            for (int j = 0; j < pts; j++) {
                Vec x(indf(i), indf(j));
                x.dimension = 2;
                x.n = n;
                points.push_back(x);
                values.push_back(func_d(x));
            }
        }
        ninds.push_back((n-1) * pts * pts);
    }
    Vec test_pnt = points[2];
    float test_val = real(values[2])(0);

    CMData data(points, values, 2, false, true, false, false);

    if (ninds != data.n_inds) return false;

    Vec data_pnt = data.points[2];
    float data_val = real(data.values[2])(0);

    if ((test_pnt - data_pnt).norm() > 1e-6) return false;
    if (abs(test_val - data_val) > 1e-6) return false;

    save(data, "testfile");
    CMData loaddata("testfile");

    Vec load_pnt = loaddata.points[2];
    float load_val = real(loaddata.values[2])(0);

    if ((test_pnt - load_pnt).norm() > 1e-6) return false;
    if (abs(test_val - load_val) > 1e-6) return false;

    return true;
}

bool readwrite_3d_with_n() {
    vector<Vec> points;
    vector<complex<Vec>> values;
    vector<int> ninds;
    for (int n = 1; n <= npts; n++) {
        for (int i = 0; i < pts; i++) {
            for (int j = 0; j < pts; j++) {
                for (int k = 0; k < pts; k++) {
                    Vec x(indf(i), indf(j), indf(k));
                    x.dimension = 2;
                    x.n = n;
                    points.push_back(x);
                    values.push_back(func_d(x));
                }
            }
        }
        ninds.push_back((n-1) * pts * pts * pts);
    }
    Vec test_pnt = points[2];
    float test_val = real(values[2])(0);

    CMData data(points, values, 3, false, true, false, false);

    if (ninds != data.n_inds) return false;

    Vec data_pnt = data.points[2];
    float data_val = real(data.values[2])(0);

    if ((test_pnt - data_pnt).norm() > 1e-6) return false;
    if (abs(test_val - data_val) > 1e-6) return false;

    save(data, "testfile");
    CMData loaddata("testfile");

    Vec load_pnt = loaddata.points[2];
    float load_val = real(loaddata.values[2])(0);

    if ((test_pnt - load_pnt).norm() > 1e-6) return false;
    if (abs(test_val - load_val) > 1e-6) return false;

    return true;
}

bool readwrite_0d_with_w_with_n() {
    vector<Vec> points;
    vector<complex<Vec>> values;
    vector<int> ninds;
    for (int n = 1; n <= npts; n++) {
        for (int i = 0; i < pts; i++) {
            Vec x(indf(i));
            x.n = n;
            points.push_back(x);
            values.push_back(func_d(x));
        }
        ninds.push_back((n-1) * pts );
    }
    Vec test_pnt = points[2];
    float test_val = real(values[2])(0);

    CMData data(points, values, 0, true, true, false, false);

    if (ninds != data.n_inds) return false;

    Vec data_pnt = data.points[2];
    float data_val = real(data.values[2])(0);

    if ((test_pnt - data_pnt).norm() > 1e-6) return false;
    if (abs(test_val - data_val) > 1e-6) return false;

    save(data, "testfile");
    CMData loaddata("testfile");

    Vec load_pnt = loaddata.points[2];
    float load_val = real(loaddata.values[2])(0);

    if ((test_pnt - load_pnt).norm() > 1e-6) return false;
    if (abs(test_val - load_val) > 1e-6) return false;

    return true;
}

bool readwrite_1d_with_w_with_n() {
    vector<Vec> points;
    vector<complex<Vec>> values;
    vector<int> ninds;
    for (int n = 1; n <= npts; n++) {
        for (int i = 0; i < pts; i++) {
            for (int j = 0; j < pts; j++) {
                Vec x(indf(i), indf(j));
                x.n = n;
                points.push_back(x);
                values.push_back(func_d(x));
            }
        }
        ninds.push_back((n-1) * pts * pts);
    }
    Vec test_pnt = points[2];
    float test_val = real(values[2])(0);

    CMData data(points, values, 1, true, true, false, false);

    if (ninds != data.n_inds) return false;

    Vec data_pnt = data.points[2];
    float data_val = real(data.values[2])(0);

    if ((test_pnt - data_pnt).norm() > 1e-6) return false;
    if (abs(test_val - data_val) > 1e-6) return false;

    save(data, "testfile");
    CMData loaddata("testfile");

    Vec load_pnt = loaddata.points[2];
    float load_val = real(loaddata.values[2])(0);

    if ((test_pnt - load_pnt).norm() > 1e-6) return false;
    if (abs(test_val - load_val) > 1e-6) return false;

    return true;
}

bool readwrite_2d_with_w_with_n() {
    vector<Vec> points;
    vector<complex<Vec>> values;
    vector<int> ninds;
    for (int n = 1; n <= npts; n++) {
        for (int i = 0; i < pts; i++) {
            for (int j = 0; j < pts; j++) {
                for (int k = 0; k < pts; k++) {
                    Vec x(indf(i), indf(j), indf(k));
                    x.n = n;
                    points.push_back(x);
                    values.push_back(func_d(x));
                }
            }
        }
        ninds.push_back((n-1) * pts * pts * pts);
    }
    Vec test_pnt = points[2];
    float test_val = real(values[2])(0);

    CMData data(points, values, 2, true, true, false, false);

    if (ninds != data.n_inds) return false;

    Vec data_pnt = data.points[2];
    float data_val = real(data.values[2])(0);

    if ((test_pnt - data_pnt).norm() > 1e-6) return false;
    if (abs(test_val - data_val) > 1e-6) return false;

    save(data, "testfile");
    CMData loaddata("testfile");

    Vec load_pnt = loaddata.points[2];
    float load_val = real(loaddata.values[2])(0);

    if ((test_pnt - load_pnt).norm() > 1e-6) return false;
    if (abs(test_val - load_val) > 1e-6) return false;

    return true;
}

bool readwrite_3d_with_w_with_n() {
    vector<Vec> points;
    vector<complex<Vec>> values;
    vector<int> ninds;
    for (int n = 1; n <= npts; n++) {
        for (int i = 0; i < pts; i++) {
            for (int j = 0; j < pts; j++) {
                for (int k = 0; k < pts; k++) {
                    for (int l = 0; l < pts; l++) {
                        Vec x(indf(i), indf(j), indf(k), indf(l));
                        x.n = n;
                        points.push_back(x);
                        values.push_back(func_d(x));
                    }
                }
            }
        }
        ninds.push_back((n-1) * pts * pts * pts * pts);
    }
    Vec test_pnt = points[2];
    float test_val = real(values[2])(0);

    CMData data(points, values, 3, true, true, false, false);

    if (ninds != data.n_inds) return false;

    Vec data_pnt = data.points[2];
    float data_val = real(data.values[2])(0);

    if ((test_pnt - data_pnt).norm() > 1e-6) return false;
    if (abs(test_val - data_val) > 1e-6) return false;

    save(data, "testfile");
    CMData loaddata("testfile");

    Vec load_pnt = loaddata.points[2];
    float load_val = real(loaddata.values[2])(0);

    if ((test_pnt - load_pnt).norm() > 1e-6) return false;
    if (abs(test_val - load_val) > 1e-6) return false;

    return true;
}
