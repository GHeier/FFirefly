#include "../cmf.hpp"
#include "../vec.hpp"

#include "../../config/load/c_config.h"

using namespace std;

int pnts = 11;

complex<Vec> func(Vec point) {
    return complex<Vec>(point.norm()*point.norm());
}

bool test_2d() {
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

    bool frompoint = test1 && test2 && test3;
    save_CMF_to_file("temp.dat", field);

    CMF filefield = load_CMF_from_file("temp.dat");

    float r4 = filefield(Vec(0.5, 0.5)).real()(0);
    float r5 = filefield(Vec(1.0, 1.0)).real()(0);
    float r6 = filefield(Vec(0.0, 0.0)).real()(0);

    bool fromfile = fabs(r4 - 0.5) < 1e-2 && fabs(r5 - 2.0) < 1e-2 && fabs(r6 - 0.0) < 1e-2;

    remove("temp.dat");
    return frompoint && fromfile;
}

bool test_2d_complex() {
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

    save_CMF_to_file("temp.dat", field);

    CMF filefield = load_CMF_from_file("temp.dat");

    Vec rv4 = filefield(Vec(0.5, 0.5)).real();
    Vec iv4 = filefield(Vec(0.5, 0.5)).imag();
    Vec rv5 = filefield(Vec(1.0, 1.0)).real();
    Vec iv5 = filefield(Vec(1.0, 1.0)).imag();
    Vec rv6 = filefield(Vec(0.0, 0.0)).real();
    Vec iv6 = filefield(Vec(0.0, 0.0)).imag();

    complex<float> r4 = complex<float>(rv4(0), iv4(0));
    complex<float> r5 = complex<float>(rv5(0), iv5(0));
    complex<float> r6 = complex<float>(rv6(0), iv6(0));

    bool test4 = fabs(r4.real() - 0.5) < 1e-2;
    bool test5 = fabs(r5.real() - 2.0) < 1e-2;
    bool test6 = fabs(r6.real() - 0.0) < 1e-2;

    remove("temp.dat");
    return test1 && test2 && test3 && test4 && test5 && test6;
}

bool test_3d() {
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

    save_CMF_to_file("temp.dat", field);

    CMF filefield = load_CMF_from_file("temp.dat");

    float r4 = filefield(Vec(0.5, 0.5, 0.5)).real()(0);
    float r5 = filefield(Vec(1.0, 1.0, 1.0)).real()(0);
    float r6 = filefield(Vec(0.0, 0.0, 0.0)).real()(0);

    bool test4 = fabs(r4 - 0.75) < 1e-2;
    bool test5 = fabs(r5 - 3.0) < 1e-2;
    bool test6 = fabs(r6 - 0.0) < 1e-2;

    remove("temp.dat");
    return test1 && test2 && test3 && test4 && test5 && test6;
}

bool test_3d_complex() {
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

    save_CMF_to_file("temp.dat", field);

    CMF filefield = load_CMF_from_file("temp.dat");

    Vec rv4 = filefield(Vec(0.5, 0.5, 0.5)).real();
    Vec iv4 = filefield(Vec(0.5, 0.5, 0.5)).imag();
    Vec rv5 = filefield(Vec(1.0, 1.0, 1.0)).real();
    Vec iv5 = filefield(Vec(1.0, 1.0, 1.0)).imag();
    Vec rv6 = filefield(Vec(0.0, 0.0, 0.0)).real();
    Vec iv6 = filefield(Vec(0.0, 0.0, 0.0)).imag();

    complex<float> r4 = complex<float>(rv4(0), iv4(0));
    complex<float> r5 = complex<float>(rv5(0), iv5(0));
    complex<float> r6 = complex<float>(rv6(0), iv6(0));

    bool test4 = fabs(r4.real() - 0.75) < 1e-2;
    bool test5 = fabs(r5.real() - 3.0) < 1e-2;
    bool test6 = fabs(r6.real() - 0.0) < 1e-2;

    remove("temp.dat");
    return test1 && test2 && test3 && test4 && test5 && test6;
}

bool test_2d_with_w() {
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

    save_CMF_to_file("temp.dat", field);

    CMF filefield = load_CMF_from_file("temp.dat");

    float r4 = filefield(Vec(0.5, 0.5), 0.5).real()(0);
    float r5 = filefield(Vec(1.0, 1.0), 1.0).real()(0);
    float r6 = filefield(Vec(0.0, 0.0), 0.0).real()(0);

    bool test4 = fabs(r4 - 0.5) < 1e-2;
    bool test5 = fabs(r5 - 2.0) < 1e-2;
    bool test6 = fabs(r6 - 0.0) < 1e-2;

    remove("temp.dat");
    return test1 && test2 && test3 && test4 && test5 && test6;
}

bool test_3d_with_w() {
    vector<Vec> points;
    vector<complex<Vec>> values;
    int dimension = 3;
    bool is_complex = false;
    bool is_vector = false;
    bool with_w = true;

    for (int i = 0; i < pnts; i++) {
        for (int j = 0; j < pnts; j++) {
            for (int k = 0; k < pnts; k++) {
                float x = 1.0 * i / (pnts - 1);
                float y = 1.0 * j / (pnts - 1);
                float w = 1.0 * k / (pnts - 1);
                Vec point(x, y, w);
                point.dimension = dimension;
                points.push_back(point);
                values.push_back(func(point));
            }
        }
    }
    auto field = CMF(points, values, dimension, with_w, is_complex, is_vector);
    float r1 = field(Vec(0.5, 0.5, 0.5), 0.5).real()(0);
    float r2 = field(Vec(1.0, 1.0, 1.0), 1.0).real()(0);
    float r3 = field(Vec(0.0, 0.0, 0.0), 0.0).real()(0);

    bool test1 = fabs(r1 - 0.75) < 1e-2;
    bool test2 = fabs(r2 - 3.0) < 1e-2;
    bool test3 = fabs(r3 - 0.0) < 1e-2;

    save_CMF_to_file("temp.dat", field);

    CMF filefield = load_CMF_from_file("temp.dat");

    float r4 = filefield(Vec(0.5, 0.5, 0.5), 0.5).real()(0);
    float r5 = filefield(Vec(1.0, 1.0, 1.0), 1.0).real()(0);
    float r6 = filefield(Vec(0.0, 0.0, 0.0), 0.0).real()(0);

    bool test4 = fabs(r4 - 0.75) < 1e-2;
    bool test5 = fabs(r5 - 3.0) < 1e-2;
    bool test6 = fabs(r6 - 0.0) < 1e-2;

    remove("temp.dat");
    return test1 && test2 && test3 && test4 && test5 && test6;
}

bool test_4d_with_w() {
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

    save_CMF_to_file("temp.dat", field);

    CMF filefield = load_CMF_from_file("temp.dat");

    float r4 = filefield(Vec(0.5, 0.5, 0.5), 0.5).real()(0);
    float r5 = filefield(Vec(1.0, 1.0, 1.0), 1.0).real()(0);
    float r6 = filefield(Vec(0.0, 0.0, 0.0), 0.0).real()(0);

    bool test4 = fabs(r4 - 1.0) < 1e-2;
    bool test5 = fabs(r5 - 4.0) < 1e-2;
    bool test6 = fabs(r6 - 0.0) < 1e-2;

    remove("temp.dat");
    return test1 && test2 && test3 && test4 && test5 && test6;
}

bool cmf_tests() {
    printf("Running CMF tests\n");
    int num_tests = 7;
    bool all_tests[num_tests] = {
        test_2d(),
        test_2d_complex(),
        test_3d(),
        test_3d_complex(),
        test_2d_with_w(),
        test_3d_with_w(),
        test_4d_with_w()
    };
    return print_test_results(all_tests, num_tests, "CMF tests");
}
