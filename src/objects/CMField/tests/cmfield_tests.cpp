#include "../../vec.hpp"
#include "../cmfield.hpp"
#include "../fields.hpp"

#include "../../../config/load/c_config.h"

using namespace std;

int pnts = 3;

complex<Vec> func(Vec point) {
    Vec p(point.norm() * point.norm());
    return complex<Vec>(p, Vec());
}

complex<Vec> func_linear(Vec point) {
    Vec p(point(0) + point(1) + point(2) + point(3));
    return complex<Vec>(p, Vec());
}

bool interp_test_2d() {
    vector<Vec> points;
    vector<complex<Vec>> values;
    int dimension = 2;
    bool is_complex = false;
    bool is_vector = false;
    bool with_w = false;
    bool with_n = false;
    printf("...\n");

    for (int i = 0; i <= pnts; i++) {
        for (int j = 0; j <= pnts; j++) {
            float x = 2.0 * (i - pnts / 2.0) / (pnts - 1);
            float y = 2.0 * (j - pnts / 2.0) / (pnts - 1);
            Vec point(x, y);
            point.dimension = dimension;
            points.push_back(point);
            values.push_back(func_linear(point));
        }
    }
    auto field = CMField(points, values, dimension, with_w, with_n, is_complex,
                         is_vector);
    float r1 = field(Vec(0.5, 0.5)).real()(0);
    float r2 = field(Vec(1.0, 1.0)).real()(0);
    float r3 = field(Vec(0.0, 0.0)).real()(0);

    bool test1 = fabs(r1 - 1.0) < 1e-2;
    bool test2 = fabs(r2 - 2.0) < 1e-2;
    bool test3 = fabs(r3 - 0.0) < 1e-2;

    printf("\nExpected: 1.0, 2.0, 0.0\n");
    printf("Got: %f, %f, %f\n\n", r1, r2, r3);

    return test1 && test2 && test3;
}

bool interp_test_2d_complex() {
    // printf("Running 2d_complex()\n");
    vector<Vec> points;
    vector<complex<Vec>> values;
    int dimension = 2;
    bool is_complex = true;
    bool is_vector = false;
    bool with_w = false;
    bool with_n = false;

    for (int i = 0; i <= pnts; i++) {
        for (int j = 0; j <= pnts; j++) {
            float x = 2.0 * (i - pnts / 2.0) / (pnts - 1);
            float y = 2.0 * (j - pnts / 2.0) / (pnts - 1);
            Vec point(x, y);
            point.dimension = dimension;
            points.push_back(point);
            values.push_back(func_linear(point));
        }
    }
    auto field = CMField(points, values, dimension, with_w, with_n, is_complex,
                         is_vector);
    Vec rv1 = field(Vec(0.5, 0.5)).real();
    Vec iv1 = field(Vec(0.5, 0.5)).imag();
    Vec rv2 = field(Vec(1.0, 1.0)).real();
    Vec iv2 = field(Vec(1.0, 1.0)).imag();
    Vec rv3 = field(Vec(0.0, 0.0)).real();
    Vec iv3 = field(Vec(0.0, 0.0)).imag();

    complex<float> r1 = complex<float>(rv1(0), iv1(0));
    complex<float> r2 = complex<float>(rv2(0), iv2(0));
    complex<float> r3 = complex<float>(rv3(0), iv3(0));

    bool test1 = fabs(r1.real() - 1.0) < 1e-2;
    bool test2 = fabs(r2.real() - 2.0) < 1e-2;
    bool test3 = fabs(r3.real() - 0.0) < 1e-2;

    printf("Expected: 1.0, 2.0, 0.0\n");
    printf("Got: %f, %f, %f\n\n", r1.real(), r2.real(), r3.real());

    return test1 && test2 && test3;
}

bool interp_test_3d() {
    // printf("Running 3d()\n");
    vector<Vec> points;
    vector<complex<Vec>> values;
    int dimension = 3;
    bool is_complex = false;
    bool is_vector = false;
    bool with_w = false;
    bool with_n = false;

    for (int i = 0; i <= pnts; i++) {
        for (int j = 0; j <= pnts; j++) {
            for (int k = 0; k <= pnts; k++) {
                float x = 2.0 * (i - pnts / 2.0) / (pnts - 1);
                float y = 2.0 * (j - pnts / 2.0) / (pnts - 1);
                float z = 2.0 * (k - pnts / 2.0) / (pnts - 1);
                Vec point(x, y, z);
                point.dimension = dimension;
                points.push_back(point);
                values.push_back(func_linear(point));
            }
        }
    }
    auto field = CMField(points, values, dimension, with_w, with_n, is_complex,
                         is_vector);
    float r1 = field(Vec(0.5, 0.5, 0.5)).real()(0);
    float r2 = field(Vec(1.0, 1.0, 1.0)).real()(0);
    float r3 = field(Vec(0.0, 0.0, 0.0)).real()(0);

    printf("\nExpected: 1.5, 3.0, 0.0\n");
    printf("Got: %f, %f, %f\n\n", r1, r2, r3);

    bool test1 = fabs(r1 - 1.5) < 1e-2;
    bool test2 = fabs(r2 - 3.0) < 1e-2;
    bool test3 = fabs(r3 - 0.0) < 1e-2;

    return test1 && test2 && test3;
}

bool interp_test_3d_complex() {
    vector<Vec> points;
    vector<complex<Vec>> values;
    int dimension = 3;
    bool is_complex = true;
    bool is_vector = false;
    bool with_w = false;
    bool with_n = false;

    for (int i = 0; i <= pnts; i++) {
        for (int j = 0; j <= pnts; j++) {
            for (int k = 0; k <= pnts; k++) {
                float x = 2.0 * (i - pnts / 2.0) / (pnts - 1);
                float y = 2.0 * (j - pnts / 2.0) / (pnts - 1);
                float z = 2.0 * (k - pnts / 2.0) / (pnts - 1);
                Vec point(x, y, z);
                point.dimension = dimension;
                points.push_back(point);
                values.push_back(func_linear(point));
            }
        }
    }
    auto field = CMField(points, values, dimension, with_w, with_n, is_complex,
                         is_vector);
    Vec rv1 = field(Vec(0.5, 0.5, 0.5)).real();
    Vec iv1 = field(Vec(0.5, 0.5, 0.5)).imag();
    Vec rv2 = field(Vec(1.0, 1.0, 1.0)).real();
    Vec iv2 = field(Vec(1.0, 1.0, 1.0)).imag();
    Vec rv3 = field(Vec(0.0, 0.0, 0.0)).real();
    Vec iv3 = field(Vec(0.0, 0.0, 0.0)).imag();

    complex<float> r1 = complex<float>(rv1(0), iv1(0));
    complex<float> r2 = complex<float>(rv2(0), iv2(0));
    complex<float> r3 = complex<float>(rv3(0), iv3(0));

    printf("Expected: 1.5, 3.0, 0.0\n");

    bool test1 = fabs(r1.real() - 1.5) < 1e-2;
    bool test2 = fabs(r2.real() - 3.0) < 1e-2;
    bool test3 = fabs(r3.real() - 0.0) < 1e-2;
    printf("Got: %f, %f, %f\n", r1.real(), r2.real(), r3.real());

    return test1 && test2 && test3;
}

bool interp_test_1d_with_w() {
    vector<Vec> points;
    vector<complex<Vec>> values;
    int dimension = 1;
    bool is_complex = false;
    bool is_vector = false;
    bool with_w = true;
    bool with_n = false;

    for (int i = 0; i <= pnts; i++) {
        for (int j = 0; j <= pnts; j++) {
            float x = 2.0 * (i - pnts / 2.0) / (pnts - 1);
            float w = 2.0 * (j - pnts / 2.0) / (pnts - 1);
            Vec point(x, w);
            point.dimension = dimension + 1;
            points.push_back(point);
            values.push_back(func_linear(point));
        }
    }
    auto field = CMField(points, values, dimension, with_w, with_n, is_complex,
                         is_vector);
    float r1 = field(Vec(0.5), 0.5).real()(0);
    float r2 = field(Vec(1.0), 1.0).real()(0);
    float r3 = field(Vec(0.0), 0.0).real()(0);

    // printf("\nExpected: 1.0, 2.0, 0.0\n");
    // printf("Got: %f, %f, %f\n\n", r1, r2, r3);

    bool test1 = fabs(r1 - 1.0) < 1e-2;
    bool test2 = fabs(r2 - 2.0) < 1e-2;
    bool test3 = fabs(r3 - 0.0) < 1e-2;

    return test1 && test2 && test3;
}

bool interp_test_2d_with_w() {
    vector<Vec> points;
    vector<complex<Vec>> values;
    int dimension = 2;
    bool is_complex = false;
    bool is_vector = false;
    bool with_w = true;
    bool with_n = false;

    for (int i = 0; i <= pnts; i++) {
        for (int j = 0; j <= pnts; j++) {
            for (int k = 0; k <= pnts; k++) {
                float x = 2.0 * (i - pnts / 2.0) / (pnts - 1);
                float y = 2.0 * (j - pnts / 2.0) / (pnts - 1);
                float w = 2.0 * (k - pnts / 2.0) / (pnts - 1);
                Vec point(x, y, w);
                point.dimension = dimension + 1;
                points.push_back(point);
                values.push_back(func_linear(point));
            }
        }
    }
    auto field = CMField(points, values, dimension, with_w, with_n, is_complex,
                         is_vector);
    float r1 = field(Vec(0.5, 0.5), 0.5).real()(0);
    float r2 = field(Vec(1.0, 1.0), 1.0).real()(0);
    float r3 = field(Vec(0.0, 0.0), 0.0).real()(0);

    // printf("\nExpected: 1.5, 3.0, 0.0\n");
    // printf("Got: %f, %f, %f\n\n", r1, r2, r3);

    bool test1 = fabs(r1 - 1.5) < 1e-2;
    bool test2 = fabs(r2 - 3.0) < 1e-2;
    bool test3 = fabs(r3 - 0.0) < 1e-2;

    return test1 && test2 && test3;
}

bool interp_test_3d_with_w() {
    vector<Vec> points;
    vector<complex<Vec>> values;
    int dimension = 3;
    bool is_complex = false;
    bool is_vector = false;
    bool with_w = true;
    bool with_n = false;

    for (int i = 0; i <= pnts; i++) {
        for (int j = 0; j <= pnts; j++) {
            for (int k = 0; k <= pnts; k++) {
                for (int l = 0; l <= pnts; l++) {
                    float x = 2.0 * (i - pnts / 2.0) / (pnts - 1);
                    float y = 2.0 * (j - pnts / 2.0) / (pnts - 1);
                    float z = 2.0 * (k - pnts / 2.0) / (pnts - 1);
                    float w = 2.0 * (l - pnts / 2.0) / (pnts - 1);
                    Vec point(x, y, z, w);
                    point.dimension = dimension + 1;
                    points.push_back(point);
                    values.push_back(func_linear(point));
                }
            }
        }
    }
    auto field = CMField(points, values, dimension, with_w, with_n, is_complex,
                         is_vector);
    float r1 = field(Vec(0.5, 0.5, 0.5), 0.5).real()(0);
    float r2 = field(Vec(1.0, 1.0, 1.0), 1.0).real()(0);
    float r3 = field(Vec(0.0, 0.0, 0.0), 0.0).real()(0);

    // printf("\nExpected: 1.0, 4.0, 0.0\n");
    // printf("Got: %f, %f, %f\n\n", r1, r2, r3);

    bool test1 = fabs(r1 - 2.0) < 1e-2;
    bool test2 = fabs(r2 - 4.0) < 1e-2;
    bool test3 = fabs(r3 - 0.0) < 1e-2;

    return test1 && test2 && test3;
}

bool interp_test_0d_with_w() {
    vector<Vec> points;
    vector<complex<Vec>> values;
    int dimension = 0;
    bool is_complex = false;
    bool is_vector = false;
    bool with_w = true;
    bool with_n = false;

    for (int i = 0; i <= pnts; i++) {
        float w = 2.0 * (i - pnts / 2.0) / (pnts - 1);
        Vec point(w);
        point.dimension = dimension;
        points.push_back(point);
        values.push_back(func_linear(point));
    }
    auto field = CMField(points, values, dimension, with_w, with_n, is_complex,
                         is_vector);
    float r1 = field(0.5).real()(0);
    float r2 = field(1.0).real()(0);
    float r3 = field(0.0).real()(0);

    // printf("\nExpected: 0.5, 1.0, 0.0\n");
    // printf("Got: %f, %f, %f\n\n", r1, r2, r3);

    bool test1 = fabs(r1 - 0.5) < 1e-2;
    bool test2 = fabs(r2 - 1.0) < 1e-2;
    bool test3 = fabs(r3 - 0.0) < 1e-2;

    return test1 && test2 && test3;
}

bool interp_test_1d_with_n() {
    vector<Vec> points;
    vector<complex<Vec>> values;
    int dimension = 1;
    bool is_complex = false;
    bool is_vector = false;
    bool with_w = false;
    bool with_n = true;

    for (int i = 1; i <= pnts; i++) {
        for (int j = 0; j <= pnts; j++) {
            float x = 2.0 * (j - pnts / 2.0) / (pnts - 1);
            Vec point(x);
            point.dimension = dimension;
            point.n = i;
            points.push_back(point);
            values.push_back(func_linear(point));
        }
    }
    auto field = CMField(points, values, dimension, with_w, with_n, is_complex,
                         is_vector);
    float r1 = field(1, Vec(0.5)).real()(0);
    float r2 = field(1, Vec(1.0)).real()(0);
    float r3 = field(1, Vec(0.0)).real()(0);

    // printf("\nExpected: 0.5, 1.0, 0.0\n");
    // printf("Got: %f, %f, %f\n\n", r1, r2, r3);

    bool test1 = fabs(r1 - 0.5) < 1e-2;
    bool test2 = fabs(r2 - 1.0) < 1e-2;
    bool test3 = fabs(r3 - 0.0) < 1e-2;

    return test1 && test2 && test3;
}

bool interp_test_1d_vector() {
    vector<Vec> points;
    vector<complex<Vec>> values;
    int dimension = 1;
    bool is_complex = false;
    bool is_vector = true;
    bool with_w = false;
    bool with_n = true;

    for (int i = 1; i <= pnts; i++) {
        for (int j = 0; j <= pnts; j++) {
            float x = 2.0 * (j - pnts / 2.0) / (pnts - 1);
            Vec point(x);
            point.dimension = dimension;
            point.n = i;
            points.push_back(point);
            values.push_back(func_linear(point));
        }
    }
    auto field = CMField(points, values, dimension, with_w, with_n, is_complex,
                         is_vector);
    Vec r1 = field(1, Vec(0.5)).real();
    Vec r2 = field(1, Vec(1.0)).real();
    Vec r3 = field(1, Vec(0.0)).real();

    // printf("\nExpected: 0.5, 1.0, 0.0\n");
    // printf("Got: %f, %f, %f\n\n", r1, r2, r3);

    bool test1 = (r1 - Vec(0.5)).norm() < 1e-2;
    bool test2 = (r2 - Vec(1.0)).norm() < 1e-2;
    bool test3 = (r3 - Vec(0.0)).norm() < 1e-2;

    return test1 && test2 && test3;
}

bool interp_test_2d_vector() {
    vector<Vec> points;
    vector<complex<Vec>> values;
    int dimension = 2;
    bool is_complex = false;
    bool is_vector = true;
    bool with_w = false;
    bool with_n = true;

    for (int i = 1; i <= pnts; i++) {
        for (int j = 0; j <= pnts; j++) {
            for (int k = 0; k <= pnts; k++) {
                float x = 2.0 * (j - pnts / 2.0) / (pnts - 1);
                float y = 2.0 * (k - pnts / 2.0) / (pnts - 1);
                Vec point(x, y);
                point.dimension = dimension;
                point.n = i;
                points.push_back(point);
                values.push_back(func_linear(point));
            }
        }
    }
    auto field = CMField(points, values, dimension, with_w, with_n, is_complex,
                         is_vector);
    Vec r1 = field(1, Vec(0.5, 0.5)).real();
    Vec r2 = field(1, Vec(1.0, 1.0)).real();
    Vec r3 = field(1, Vec(0.0, 0.0)).real();

    // printf("\nExpected: 0.5, 1.0, 0.0\n");
    // printf("Got: %f, %f, %f\n\n", r1, r2, r3);

    bool test1 = (r1 - Vec(1.0, 0.0)).norm() < 1e-2;
    bool test2 = (r2 - Vec(2.0, 0.0)).norm() < 1e-2;
    bool test3 = (r3 - Vec(0.0, 0.0)).norm() < 1e-2;

    return test1 && test2 && test3;
}

bool interp_test_3d_vector() {
    vector<Vec> points;
    vector<complex<Vec>> values;
    int dimension = 3;
    bool is_complex = false;
    bool is_vector = true;
    bool with_w = false;
    bool with_n = true;

    for (int i = 1; i <= pnts; i++) {
        for (int j = 0; j <= pnts; j++) {
            for (int k = 0; k <= pnts; k++) {
                for (int l = 0; l <= pnts; l++) {
                    float x = 2.0 * (j - pnts / 2.0) / (pnts - 1);
                    float y = 2.0 * (k - pnts / 2.0) / (pnts - 1);
                    float z = 2.0 * (l - pnts / 2.0) / (pnts - 1);
                    Vec point(x, y, z);
                    point.dimension = dimension;
                    point.n = i;
                    points.push_back(point);
                    values.push_back(func_linear(point));
                }
            }
        }
    }
    auto field = CMField(points, values, dimension, with_w, with_n, is_complex,
                         is_vector);
    Vec r1 = field(1, Vec(0.5, 0.5, 0.5)).real();
    Vec r2 = field(1, Vec(1.0, 1.0, 1.0)).real();
    Vec r3 = field(1, Vec(0.0, 0.0, 0.0)).real();

    // printf("\nExpected: 0.5, 1.0, 0.0\n");
    // printf("Got: %f, %f, %f\n\n", r1, r2, r3);
    // cout << "Got:\n" << r1 << "\n" << r2 << "\n" << r3 << endl;

    bool test1 = (r1 - Vec(1.5, 0.0)).norm() < 1e-2;
    bool test2 = (r2 - Vec(3.0, 0.0)).norm() < 1e-2;
    bool test3 = (r3 - Vec(0.0, 0.0)).norm() < 1e-2;

    return test1 && test2 && test3;
}

bool interp_test_1d() {
    vector<Vec> points;
    vector<complex<Vec>> values;
    int dimension = 1;
    bool is_complex = false;
    bool is_vector = false;
    bool with_w = false;
    bool with_n = false;

    for (int i = 1; i <= pnts; i++) {
        float x = 2.0 * (i - pnts / 2.0) / (pnts - 1);
        Vec point(x);
        point.dimension = dimension;
        point.n = i;
        points.push_back(point);
        values.push_back(func_linear(point));
    }
    auto field = CMField(points, values, dimension, with_w, with_n, is_complex,
                         is_vector);
    Vec r1 = field(Vec(0.5)).real();
    Vec r2 = field(Vec(1.0)).real();
    Vec r3 = field(Vec(0.0)).real();

    // printf("\nExpected: 0.5, 1.0, 0.0\n");
    // printf("Got: %f, %f, %f\n\n", r1, r2, r3);
    // cout << "Got:\n" << r1 << "\n" << r2 << "\n" << r3 << endl;

    bool test1 = (r1 - Vec(0.5, 0.0)).norm() < 1e-2;
    bool test2 = (r2 - Vec(1.0, 0.0)).norm() < 1e-2;
    bool test3 = (r3 - Vec(0.0, 0.0)).norm() < 1e-2;

    return test1 && test2 && test3;
}

bool interp_test_1d_complex() {
    vector<Vec> points;
    vector<complex<Vec>> values;
    int dimension = 1;
    bool is_complex = true;
    bool is_vector = false;
    bool with_w = false;
    bool with_n = false;

    for (int i = 1; i <= pnts; i++) {
        float x = 2.0 * (i - pnts / 2.0) / (pnts - 1);
        Vec point(x);
        point.dimension = dimension;
        point.n = i;
        points.push_back(point);
        values.push_back(func_linear(point));
    }
    auto field = CMField(points, values, dimension, with_w, with_n, is_complex,
                         is_vector);
    Vec r1 = field(Vec(0.5)).real();
    Vec r2 = field(Vec(1.0)).real();
    Vec r3 = field(Vec(0.0)).real();

    // printf("\nExpected: 0.5, 1.0, 0.0\n");
    // printf("Got: %f, %f, %f\n\n", r1, r2, r3);
    // cout << "Got:\n" << r1 << "\n" << r2 << "\n" << r3 << endl;

    bool test1 = (r1 - Vec(0.5, 0.0)).norm() < 1e-2;
    bool test2 = (r2 - Vec(1.0, 0.0)).norm() < 1e-2;
    bool test3 = (r3 - Vec(0.0, 0.0)).norm() < 1e-2;

    return test1 && test2 && test3;
}

bool cmfield_tests() {
    int num_tests = 14;
    bool all_tests[num_tests] = {
        interp_test_2d(),        interp_test_2d_complex(),
        interp_test_3d(),        interp_test_3d_complex(),
        interp_test_0d_with_w(), interp_test_1d_with_w(),
        interp_test_2d_with_w(), interp_test_3d_with_w(),
        interp_test_1d_with_n(), interp_test_1d_vector(),
        interp_test_2d_vector(), interp_test_3d_vector(),
        interp_test_1d(),        interp_test_1d_complex(),
    };
    return print_test_results(all_tests, num_tests, "CMField tests");
}
