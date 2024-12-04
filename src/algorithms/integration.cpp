#include <vector>
#include <functional>

#include "../objects/vec.h"
#include "../objects/surfaces.h"
#include "../config/load/cpp_config.h"

int s_div = (dim == 3) ? 40 : 300; // Number of integral surface divisions
int s_pts = (dim == 3) ? 50 : 1000; // Number of integral surfaces
// Gaussian integration constants
float weights_0th[1] = {2.0}; float * w0 = weights_0th;
float weights_1st[2] = {1.0, 1.0}; float * w1 = weights_1st;
float weights_2nd[3] = {5.0/9.0, 8.0/9.0, 5.0/9.0}; float * w2 = weights_2nd;
float weights_3rd[4] = {0.347855, 0.652145, 0.652145, 0.347855}; float * w3 = weights_3rd;
float weights_4th[5] = {0.236927, 0.478629, 0.568889, 0.478629, 0.236927}; float * w4 = weights_4th;
float *weights[5] = {w0, w1, w2, w3, w4};

float points_0th[1] = {0}; float *p0 = points_0th;
float points_1st[2] = {-1/pow(3,0.5), 1/pow(3,0.5)}; float *p1 = points_1st;
float points_2nd[3] = {-pow(3/5,0.5), 0, pow(3/5,0.5)}; float *p2 = points_2nd;
float points_3rd[4] = {-0.861136, -0.339981, 0.339981, 0.861136}; float *p3 = points_3rd;
float points_4th[5] = {-0.90618, -0.538469, 0, 0.538469, 0.90618}; float *p4 = points_4th;

float *points[5] = {p0, p1, p2, p3, p4};
// Defines the constants of the integral spacing based around the upper and lower energy bounds
// found in get_bounds
void get_spacing_curve_consts(float w, float a, float b, float &A, float &upr, float &lwr) {
    A = b - w;
    lwr = (a - w) / A;
    upr = 1;
    if ( w - a > b - w) {
        A = w - a;
        lwr = -1;
        upr = (b - w) / A;
    }
}

// Creates the array of energies to be integrated over
void get_spacing_vec(vector<float> &spacing, float w, float a, float b, int pts) {
    float A, upr, lwr;
    get_spacing_curve_consts(w, a, b, A, upr, lwr);

    auto spacing_curve = [A, w] (float i, float pts) { 
        float x = -1 + 2 * i / pts;
        return A * x + w;
    };

    float r = spacing_curve(0, pts);
    for (int i = 0; r < b; i++) {
        float t = i;
        r = spacing_curve(i, pts);
        spacing.push_back(r);
    }
}

void get_surface_transformed_bounds(float &upper, float &lower, function<float(Vec)> func) {
    auto get_k = [] (float i, int pts) { return k_max*(2.0*i/(pts-1.0)-1.0); };

    upper = 0; lower = 1000;
    int pts = 100;
    for (float i = 0; i < pts; i++) {
        float x = get_k(i, pts);
        for (float j = 0; j < pts; j++) {
            float y = get_k(j, pts);
            for (float k = 0; k < pts * (dim%2) + 1 * ((dim+1)%2); k++) {
                float z = get_k(k, pts);
                Vec k_val(x, y, z);
                float val = func(k_val);
                if (val > upper) upper = val;
                if (val < lower) lower = val;
            }
        }
    }
    lower *= 0.99; upper *= 0.99;
}

float surface_transform_integral(function<float(Vec)> integrand,
        function<float(Vec)> func,
        function<float(Vec)> func_diff,
        vector<float> &svals) {
    vector<vector<float>> tetrahedrons {
        {1, 2, 3, 5}, 
        {1, 3, 4, 5},
        {2, 5, 6, 3},
        {4, 5, 8, 3},
        {5, 8, 7, 3},
        {5, 6, 7, 3}
    };

    double sum = 0;
    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < s_div; i++) {
        for (int j = 0; j < s_div; j++) {
            for (int k = 0; k < s_div * (dim%2) + 1 * ((dim+1)%2); k++) {
                vector<Vec> points = points_from_indices(func, i, j, k, s_div);
                float min = 1000, max = -1000;
                for (Vec p : points) {
                    if (func(p) < min) min = func(p);
                    if (func(p) > max) max = func(p);
                }
                pair<int, int> index_and_length = get_index_and_length(min, max, svals);
                if (index_and_length.first == -1) continue;

                for (int c = 0; c < 6; c++) {

                    vector<Vec> ep_points(4);
                    for (int p = 0; p < 4; p++) {
                        ep_points[p] = points[tetrahedrons[c][p]-1];
                    }
                    for (int x = 0; x < index_and_length.second; x++) {
                        int ind = x + index_and_length.first;
                        float s_val = svals[ind];

                        if (not surface_inside_tetrahedron(s_val, ep_points)) continue;
                        vector<Vec> corner_points = points_in_tetrahedron(func, s_val, ep_points);

                        Vec average;

                        float b = 0;
                        if (corner_points[3] == average) b = 1.0;

                        for (Vec q : corner_points) {
                            average = (q + average);
                        }
                        average = average / (4-b);

                        float A = area_in_corners(corner_points);
                        if (dim == 2) A *= s_div / (2*k_max);
                        Vec k_point = average; k_point.area = A;
                        k_point.w = s_val;

                        // Surfaces being summed over. Depending on whether the surface is 
                        // at the beginning of the end of the list, we weight it differently
                        // This is trapezoidal integration
                        double dS;
                        if (ind == 0) dS = (svals[ind+1] - svals[ind]) / 2;
                        else if (ind == svals.size()-1) dS = (svals[ind] - svals[ind-1]) / 2;
                        else dS = (svals[ind+1] - svals[ind-1]) / 2;

                        double val = integrand(k_point) * k_point.area / func_diff(k_point);
                        sum += val * dS;
                    }
                }
            }
        }
    }
    return (float)sum;
}

float trapezoidal_integration(const function<double(float, float, float)> &f, float x0, float x1, float y0, float y1, float z0, float z1, int num_points) {
    double sum = 0;
    float dx = (x1 - x0) / (num_points - 1);
    float dy = (y1 - y0) / (num_points - 1);
    float dz = (z1 - z0) / (num_points - 1);
    if (dim == 2) dz = 1;
    int counter = 0;
    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < num_points; i++) {
        float x = x0 + i*dx;
        for (int j = 0; j < num_points; j++) {
            float y = y0 + j*dy;
            for (float k = 0; k < num_points * (dim%2) + 1 * ((dim+1)%2); k++) {
                float z = z0 + k*dz;
                float weight = 1.0;
                if (i == 0 or i == num_points - 1) weight /= 2.0;
                if (j == 0 or j == num_points - 1) weight /= 2.0;
                if ( (k == 0 or k == num_points - 1) and dim == 3) weight /= 2.0;

                sum += weight * f(x,y,z) * dx * dy * dz;
            }
        }
    }
    return (float)sum;
}

complex<float> complex_trapezoidal_integration(const function<complex<float>(float, float, float)> &f, float x0, float x1, float y0, float y1, float z0, float z1, int num_points) {
    double sum_real = 0;
    double sum_imag = 0;
    float dx = (x1 - x0) / (num_points - 1);
    float dy = (y1 - y0) / (num_points - 1);
    float dz = (z1 - z0) / (num_points - 1);
    float num_z_points = num_points;
    if (dim == 2) dz = 1;

    #pragma omp parallel for reduction(+:sum_real, sum_imag)
    for (int i = 0; i < num_points; i++) {
        float x = x0 + i * dx;
        for (int j = 0; j < num_points; j++) {
            float y = y0 + j * dy;
            for (float k = 0; k < num_points; k++) {
                float z = z0 + k * dz;
                float weight = 1.0;
                if (i == 0 || i == num_points - 1) weight /= 2.0;
                if (j == 0 || j == num_points - 1) weight /= 2.0;
                if ((k == 0 || k == num_points - 1) && dim == 3) weight /= 2.0;

                double f_real = f(x, y, z).real();
                double f_imag = f(x, y, z).imag();
                sum_real += weight * f_real * dx * dy * dz;
                sum_imag += weight * f_imag * dx * dy * dz;
            }
        }
    }
    return complex<float>(sum_real, sum_imag);
}
