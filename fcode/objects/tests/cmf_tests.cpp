#include "../cmf.hpp"
#include "../vec.hpp"

#include "../../config/load/c_config.h"

using namespace std;

int pnts = 11;

complex<Vec> func(Vec point) {
    return complex<Vec>(point.norm()*point.norm());
}

bool test_2d() {
    printf("Running 2D test\n");
    vector<Vec> points;
    vector<complex<Vec>> values;
    int dimension = 2;
    bool is_complex = false;
    bool is_vector = false;
    bool with_w = false;

    for (int i = 0; i < pnts; i++) {
        for (int j = 0; j < pnts; j++) {
            float x = 1.0 * i / (pnts - 1);
            float y = 1.0 * j / (pnts - 1);
            Vec point(x, y);
            point.dimension = dimension;
            points.push_back(point);
            values.push_back(func(point));
        }
    }
    auto field = CMF(points, values, dimension, with_w, is_complex, is_vector);
    float r1 = field(Vec(0.5, 0.5)).real()(0);
    float r2 = field(Vec(1.0, 1.0)).real()(0);
    float r3 = field(Vec(0.0, 0.0)).real()(0);

    bool test1 = fabs(r1 - 0.5) < 1e-2;
    bool test2 = fabs(r2 - 2.0) < 1e-2;
    bool test3 = fabs(r3 - 0.0) < 1e-2;

    return test1 && test2 && test3;
}

bool test_2d_complex() {
    printf("Running 2D complex test\n");
    vector<Vec> points;
    vector<complex<Vec>> values;
    int dimension = 2;
    bool is_complex = true;
    bool is_vector = false;
    bool with_w = false;

    for (int i = 0; i < pnts; i++) {
        for (int j = 0; j < pnts; j++) {
            float x = 1.0 * i / (pnts - 1);
            float y = 1.0 * j / (pnts - 1);
            Vec point(x, y);
            point.dimension = dimension;
            points.push_back(point);
            values.push_back(func(point));
        }
    }
    auto field = CMF(points, values, dimension, with_w, is_complex, is_vector);
    Vec rv1 = field(Vec(0.5, 0.5)).real();
    Vec iv1 = field(Vec(0.5, 0.5)).imag();
    Vec rv2 = field(Vec(1.0, 1.0)).real();
    Vec iv2 = field(Vec(1.0, 1.0)).imag();
    Vec rv3 = field(Vec(0.0, 0.0)).real();
    Vec iv3 = field(Vec(0.0, 0.0)).imag();

    complex<float> r1 = complex<float>(rv1(0), iv1(0));
    complex<float> r2 = complex<float>(rv2(0), iv2(0));
    complex<float> r3 = complex<float>(rv3(0), iv3(0));

    bool test1 = fabs(r1.real() - 0.5) < 1e-2;
    bool test2 = fabs(r2.real() - 2.0) < 1e-2;
    bool test3 = fabs(r3.real() - 0.0) < 1e-2;

    return test1 && test2 && test3;
}

bool test_3d() {
    printf("Running 3D test\n");
    vector<Vec> points;
    vector<complex<Vec>> values;
    int dimension = 3;
    bool is_complex = false;
    bool is_vector = false;
    bool with_w = false;

    for (int i = 0; i < pnts; i++) {
        for (int j = 0; j < pnts; j++) {
            for (int k = 0; k < pnts; k++) {
                float x = 1.0 * i / (pnts - 1);
                float y = 1.0 * j / (pnts - 1);
                float z = 1.0 * k / (pnts - 1);
                Vec point(x, y, z);
                point.dimension = dimension;
                points.push_back(point);
                values.push_back(func(point));
            }
        }
    }
    auto field = CMF(points, values, dimension, with_w, is_complex, is_vector);
    float r1 = field(Vec(0.5, 0.5, 0.5)).real()(0);
    float r2 = field(Vec(1.0, 1.0, 1.0)).real()(0);
    float r3 = field(Vec(0.0, 0.0, 0.0)).real()(0);

    bool test1 = fabs(r1 - 0.75) < 1e-2;
    bool test2 = fabs(r2 - 3.0) < 1e-2;
    bool test3 = fabs(r3 - 0.0) < 1e-2;

    return test1 && test2 && test3;
}

bool test_3d_complex() {
    printf("Running 3D complex test\n");
    vector<Vec> points;
    vector<complex<Vec>> values;
    int dimension = 3;
    bool is_complex = true;
    bool is_vector = false;
    bool with_w = false;

    for (int i = 0; i < pnts; i++) {
        for (int j = 0; j < pnts; j++) {
            for (int k = 0; k < pnts; k++) {
                float x = 1.0 * i / (pnts - 1);
                float y = 1.0 * j / (pnts - 1);
                float z = 1.0 * k / (pnts - 1);
                Vec point(x, y, z);
                point.dimension = dimension;
                points.push_back(point);
                values.push_back(func(point));
            }
        }
    }
    auto field = CMF(points, values, dimension, with_w, is_complex, is_vector);
    Vec rv1 = field(Vec(0.5, 0.5, 0.5)).real();
    Vec iv1 = field(Vec(0.5, 0.5, 0.5)).imag();
    Vec rv2 = field(Vec(1.0, 1.0, 1.0)).real();
    Vec iv2 = field(Vec(1.0, 1.0, 1.0)).imag();
    Vec rv3 = field(Vec(0.0, 0.0, 0.0)).real();
    Vec iv3 = field(Vec(0.0, 0.0, 0.0)).imag();

    complex<float> r1 = complex<float>(rv1(0), iv1(0));
    complex<float> r2 = complex<float>(rv2(0), iv2(0));
    complex<float> r3 = complex<float>(rv3(0), iv3(0));

    bool test1 = fabs(r1.real() - 0.75) < 1e-2;
    bool test2 = fabs(r2.real() - 3.0) < 1e-2;
    bool test3 = fabs(r3.real() - 0.0) < 1e-2;

    return test1 && test2 && test3;
}

bool test_2d_with_w() {
    printf("Running 2D with_w test\n");
    vector<Vec> points;
    vector<complex<Vec>> values;
    int dimension = 2;
    bool is_complex = false;
    bool is_vector = false;
    bool with_w = true;

    for (int i = 0; i < pnts; i++) {
        for (int j = 0; j < pnts; j++) {
            float x = 1.0 * i / (pnts - 1);
            float w = 1.0 * j / (pnts - 1);
            Vec point(x, w);
            point.dimension = dimension;
            points.push_back(point);
            values.push_back(func(point));
        }
    }
    auto field = CMF(points, values, dimension, with_w, is_complex, is_vector);
    float r1 = field(Vec(0.5, 0.5), 0.5).real()(0);
    float r2 = field(Vec(1.0, 1.0), 1.0).real()(0);
    float r3 = field(Vec(0.0, 0.0), 0.0).real()(0);

    bool test1 = fabs(r1 - 0.5) < 1e-2;
    bool test2 = fabs(r2 - 2.0) < 1e-2;
    bool test3 = fabs(r3 - 0.0) < 1e-2;

    printf("Test 1\n Expected: 0.5\n Got: %f\n", r1);
    printf("Test 2\n Expected: 2.0\n Got: %f\n", r2);
    printf("Test 3\n Expected: 0.0\n Got: %f\n", r3);

    return test1 && test2 && test3;
}

bool test_4d_with_w() {
    printf("Running 4D test\n");
    vector<Vec> points;
    vector<complex<Vec>> values;
    int dimension = 4;
    bool is_complex = false;
    bool is_vector = false;
    bool with_w = true;

    for (int i = 0; i < pnts; i++) {
        for (int j = 0; j < pnts; j++) {
            for (int k = 0; k < pnts; k++) {
                for (int l = 0; l < pnts; l++) {
                    float x = 1.0 * i / (pnts - 1);
                    float y = 1.0 * j / (pnts - 1);
                    float z = 1.0 * k / (pnts - 1);
                    float w = 1.0 * l / (pnts - 1);
                    Vec point(x, y, z, w);
                    point.dimension = dimension;
                    points.push_back(point);
                    values.push_back(func(point));
                }
            }
        }
    }
    auto field = CMF(points, values, dimension, with_w, is_complex, is_vector);
    float r1 = field(Vec(0.5, 0.5, 0.5), 0.5).real()(0);
    float r2 = field(Vec(1.0, 1.0, 1.0), 1.0).real()(0);
    float r3 = field(Vec(0.0, 0.0, 0.0), 0.0).real()(0);

    bool test1 = fabs(r1 - 1.0) < 1e-2;
    bool test2 = fabs(r2 - 4.0) < 1e-2;
    bool test3 = fabs(r3 - 0.0) < 1e-2;

    printf("Test 1\n Expected: 1.0\n Got: %f\n", r1);
    printf("Test 2\n Expected: 4.0\n Got: %f\n", r2);
    printf("Test 3\n Expected: 0.0\n Got: %f\n", r3);

    return test1 && test2 && test3;
}

bool cmf_tests() {
    printf("Running CMF tests\n");
    int num_tests = 6;
    bool all_tests[num_tests] = {
        test_2d(),
        test_2d_complex(),
        test_3d(),
        test_3d_complex(),
        test_2d_with_w(),
        test_4d_with_w()
    };
    return print_test_results(all_tests, num_tests, "CMF tests");
}
