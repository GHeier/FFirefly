#include "../../vec.hpp"
#include "../fields.hpp"
#include "../../../config/load/c_config.h"

int fpnts = 10;

complex<Vec> fieldfunc(Vec point) {
    Vec p(point(0) + point(1) + point(2) + point(3));
    return complex<Vec>(p, Vec());
}

bool Field_R_interp_test_2d() {
    vector<Vec> points;
    vector<complex<Vec>> values;
    int dimension = 2;
    bool is_complex = false;
    bool is_vector = false;
    bool with_w = false;
    bool with_n = false;

    for (int i = 0; i <= fpnts; i++) {
        for (int j = 0; j <= fpnts; j++) {
            float x = 2.0 * (i - fpnts / 2.0) / (fpnts - 1);
            float y = 2.0 * (j - fpnts / 2.0) / (fpnts - 1);
            Vec point(x, y);
            point.dimension = dimension;
            points.push_back(point);
            values.push_back(fieldfunc(point));
        }
    }
    CMData data(points, values, dimension, with_w, with_n, is_complex, is_vector);
    save(data, "temp.dat");

    auto real_field = Field_R("temp.dat");

    float r1 = real_field(Vec(0.5, 0.5));
    float r2 = real_field(Vec(1.0, 1.0));
    float r3 = real_field(Vec(0.0, 0.0));

    float e1 = 1.0;
    float e2 = 2.0;
    float e3 = 0.0;

    bool test1 = fabs(r1 - e1) < 1e-2;
    bool test2 = fabs(r2 - e2) < 1e-2;
    bool test3 = fabs(r3 - e3) < 1e-2;


    // printf("\nExpected: 1.0, 2.0, 0.0\n");
    // printf("Got: %f, %f, %f\n\n", r1, r2, r3);

    return test1 && test2 && test3;
}

bool Field_C_interp_test_2d() {
    vector<Vec> points;
    vector<complex<Vec>> values;
    int dimension = 2;
    bool is_complex = true;
    bool is_vector = false;
    bool with_w = false;
    bool with_n = false;

    for (int i = 0; i <= fpnts; i++) {
        for (int j = 0; j <= fpnts; j++) {
            float x = 2.0 * (i - fpnts / 2.0) / (fpnts - 1);
            float y = 2.0 * (j - fpnts / 2.0) / (fpnts - 1);
            Vec point(x, y);
            point.dimension = dimension;
            points.push_back(point);
            values.push_back(fieldfunc(point));
        }
    }
    CMData data(points, values, dimension, with_w, with_n, is_complex, is_vector);
    save(data, "temp.dat");

    auto complex_field = Field_C("temp.dat");
    complex<float> r1 = complex_field(Vec(0.5, 0.5));
    complex<float> r2 = complex_field(Vec(1.0, 1.0));
    complex<float> r3 = complex_field(Vec(0.0, 0.0));

    complex<float> e1 = complex<float>(1.0, 0);
    complex<float> e2 = complex<float>(2.0, 0);
    complex<float> e3 = complex<float>(0.0, 0);



    bool test1 = fabs(r1 - e1) < 1e-2;
    bool test2 = fabs(r2 - e2) < 1e-2;
    bool test3 = fabs(r3 - e3) < 1e-2;

    // printf("\nExpected: 1.0, 2.0, 0.0\n");
    // printf("Got: %f, %f, %f\n\n", r1, r2, r3);

    return test1 && test2 && test3;
}


bool Field_C_interp_test_2d_with_w() {
    vector<Vec> points;
    vector<complex<Vec>> values;
    int dimension = 2;
    bool is_complex = true;
    bool is_vector = false;
    bool with_w = true;
    bool with_n = false;

    for (int i = 0; i <= fpnts; i++) {
        for (int j = 0; j <= fpnts; j++) {
            for (int k = 0; k <= fpnts; k++) {
                float x = 2.0 * (i - fpnts / 2.0) / (fpnts - 1);
                float y = 2.0 * (j - fpnts / 2.0) / (fpnts - 1);
                float w = 2.0 * (k - fpnts / 2.0) / (fpnts - 1);
                Vec point(x, y, w);
                point.dimension = dimension + 1;
                //cout << "p: " << point << endl;
                points.push_back(point);
                values.push_back(fieldfunc(point));
            }
        }
    }
    CMData data(points, values, dimension, with_w, with_n, is_complex, is_vector);
    save(data, "temp.dat");

    auto complex_field = Field_C("temp.dat");
    complex<float> r1 = complex_field(Vec(0.5, 0.5));
    complex<float> r2 = complex_field(Vec(1.0, 1.0));
    complex<float> r3 = complex_field(Vec(0.0, 0.0));

    complex<float> e1 = complex<float>(1.0, 0);
    complex<float> e2 = complex<float>(2.0, 0);
    complex<float> e3 = complex<float>(0.0, 0);

    bool test1 = fabs(r1 - e1) < 1e-2;
    bool test2 = fabs(r2 - e2) < 1e-2;
    bool test3 = fabs(r3 - e3) < 1e-2;

    // printf("\nExpected: 1.0, 2.0, 0.0\n");
    // printf("Got: %f, %f, %f\n\n", r1, r2, r3);

    return test1 && test2 && test3;
}

bool field_tests() {
    int num_tests = 3;
    bool all_tests[num_tests] = {
        Field_R_interp_test_2d(),
        Field_C_interp_test_2d(),
        Field_C_interp_test_2d_with_w(),
    };
    return print_test_results(all_tests, num_tests, "Field tests");
}
