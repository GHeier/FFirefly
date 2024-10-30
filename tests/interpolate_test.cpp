#include <iostream>
#include <vector>
#include <cassert>
#include "../gap/interpolate.h"

using namespace std;

void test_interpolate_1D() {
    vector<float> f = {0.0, 1.0, 2.0, 3.0}; // Linear function f(x) = x over [0, 3]
    float x_min = 0.0;
    float x_max = 3.0;

    assert(interpolate_1D(1.5, x_min, x_max, f) == 1.5); // Middle of [1, 2]
    assert(interpolate_1D(0.0, x_min, x_max, f) == 0.0); // Boundary check
    assert(interpolate_1D(3.0, x_min, x_max, f) == 3.0); // Boundary check

    cout << "1D interpolation tests passed!" << endl;
}

void test_interpolate_2D() {
    vector<vector<float>> f = {
        {0.0, 1.0},
        {1.0, 2.0}
    }; // Function f(x, y) = x + y over a 2x2 grid

    float x_min = 0.0, x_max = 1.0;
    float y_min = 0.0, y_max = 1.0;

    assert(interpolate_2D(0.5, 0.5, x_min, x_max, y_min, y_max, f) == 1.0); // Center point
    assert(interpolate_2D(0.0, 0.0, x_min, x_max, y_min, y_max, f) == 0.0); // Bottom-left corner
    assert(interpolate_2D(1.0, 1.0, x_min, x_max, y_min, y_max, f) == 2.0); // Top-right corner

    cout << "2D interpolation tests passed!" << endl;
}

void test_interpolate_3D() {
    vector<vector<vector<float>>> f = {
        {{0.0, 1.0}, {1.0, 2.0}},
        {{1.0, 2.0}, {2.0, 3.0}}
    }; // Function f(x, y, z) = x + y + z over a 2x2x2 grid

    float x_min = 0.0, x_max = 1.0;
    float y_min = 0.0, y_max = 1.0;
    float z_min = 0.0, z_max = 1.0;

    assert(interpolate_3D(0.5, 0.5, 0.5, x_min, x_max, y_min, y_max, z_min, z_max, f) == 1.5); // Center point
    assert(interpolate_3D(0.0, 0.0, 0.0, x_min, x_max, y_min, y_max, z_min, z_max, f) == 0.0); // Bottom-front-left corner
    assert(interpolate_3D(1.0, 1.0, 1.0, x_min, x_max, y_min, y_max, z_min, z_max, f) == 3.0); // Top-back-right corner

    cout << "3D interpolation tests passed!" << endl;
}

void test_interpolate_4D() {
    vector<vector<vector<vector<float>>>> f = {
        {{{0.0, 1.0}, {1.0, 2.0}}, {{1.0, 2.0}, {2.0, 3.0}}},
        {{{1.0, 2.0}, {2.0, 3.0}}, {{2.0, 3.0}, {3.0, 4.0}}}
    }; // Function f(x, y, z, w) = x + y + z + w over a 2x2x2x2 grid

    float x_min = 0.0, x_max = 1.0;
    float y_min = 0.0, y_max = 1.0;
    float z_min = 0.0, z_max = 1.0;
    float w_min = 0.0, w_max = 1.0;

    assert(interpolate_4D(0.5, 0.5, 0.5, 0.5, x_min, x_max, y_min, y_max, z_min, z_max, w_min, w_max, f) == 2.0); // Center point
    assert(interpolate_4D(0.0, 0.0, 0.0, 0.0, x_min, x_max, y_min, y_max, z_min, z_max, w_min, w_max, f) == 0.0); // Origin
    assert(interpolate_4D(1.0, 1.0, 1.0, 1.0, x_min, x_max, y_min, y_max, z_min, z_max, w_min, w_max, f) == 4.0); // Top-back-right-upper corner

    cout << "4D interpolation tests passed!" << endl;
}

