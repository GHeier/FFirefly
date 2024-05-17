#include <iostream>
#include <unistd.h>
#include <ctime>

#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <map>
#include <boost/functional/hash.hpp>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 
//#include <lambda_lanczos/lambda_lanczos.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/quadrature/gauss.hpp>

#include "calculations.h"
#include "cfg.h"
#include "fermi_surface.h"
#include "potential.h"
#include "vec.h"
#include "utilities.h"
#include "frequency_inclusion.hpp"

using std::cout;
using std::endl;
using std::sort;
using std::vector;
using std::unordered_map;
//using lambda_lanczos::LambdaLanczos;
using namespace Eigen;

vector<vector<Vec>> freq_tetrahedron_method(double mu) {
    assert( l % 2 != 0); // N must be odd that way frequencies are evenly spaced
    vector<vector<Vec>> basis;

    double points_0th[1] = {0}; double *p0 = points_0th;
    double points_1st[2] = {-1/pow(3,0.5), 1/pow(3,0.5)}; double *p1 = points_1st;
    double points_2nd[3] = {-pow(3/5,0.5), 0, pow(3/5,0.5)}; double *p2 = points_2nd;
    double points_3rd[4] = {-0.861136, -0.339981, 0.339981, 0.861136}; double *p3 = points_3rd;
    double points_4th[5] = {-0.90618, -0.538469, 0, 0.538469, 0.90618}; double *p4 = points_4th;

    double *points[5] = {p0, p1, p2, p3, p4};

    for (int i = 0; i < l; i++) {
        double ep = w_D * points[l-1][i] + mu;
        vector<Vec> layer = tetrahedron_method(e_base_avg, Vec(0,0,0), ep);
        basis.push_back(layer);
    }
    return basis;
}

MatrixXd create_P_freq(vector<vector<Vec>> &k, double T, const unordered_map<double, vector<vector<vector<double>>>> &chi_cube2) {
    int size = 0;
    for (int i = 0; i < k.size(); i++) {
        size += k[i].size();
    }

    MatrixXd P(size, size);
    cout << "Creating P Matrix with frequency\n";
    for (int i = 0; i < k.size(); i++) {

        int ind1 = 0;
        for (int temp = 0; temp < i; temp++)
            ind1 += k[temp].size();

        #pragma omp parallel for
        for (int j = 0; j < k[i].size(); j++) {
            Vec k1 = k[i][j];
            for (int x = 0; x < k.size(); x++) {

                int ind2 = 0;
                for (int temp = 0; temp < x; temp++)
                    ind2 += k[temp].size();

                for (int y = 0; y < k[x].size(); y++) {
                    Vec k2 = k[x][y];
                    double d1 = pow(k1.area/vp(k1),0.5); 
                    double d2 = pow(k2.area/vp(k2),0.5); 
                    // f * d_epsilon
                    double fde1 = f_singlet(w_D * points[l-1][i], T) * weights[l-1][i];
                    double fde2 = f_singlet(w_D * points[l-1][x], T) * weights[l-1][x];
                    double w = w_D * (points[l-1][x] - points[l-1][i]);

                    P(ind1 + j,ind2 + y) = - d1 * d2 * pow(fde1*fde2,0.5) * V(k1, k2, w, T, chi_cube2); 
                    if (ind1 + j >= size or ind2 + y >= size)
                        cout << "Error: " << ind1 + j << " " << ind2 + y << " " << size << endl;
                }
            }
        }
    }
    cout << "P Matrix Created\n";

    //return P * 2 * w_D / (l * k_size);
    return P * w_D * (2 / pow(2*M_PI, dim)); 
}

unordered_map <double, vector<vector<vector<double>>>> chi_cube_freq(double T, double mu, double DOS) {
    vector<double> des;
    for (int i = 0; i < l; i++) {
        double p1 = w_D * points[l-1][i];
        for (int j = 0; j < l; j++) {
            double p2 = w_D * points[l-1][j];
            double w = p1 - p2;
            w = round(w, 6);
            if (count(des.begin(), des.end(), w) == 0)
                des.push_back(w);
        }
    }
    unordered_map <double, vector<vector<vector<double>>>> cube_freq_map;
    for (double w : des) {
        auto cube = chi_cube(T, mu, DOS, w);
        cube_freq_map.insert(pair<double, vector<vector<vector<double>>>>(w, cube));
    }
    return cube_freq_map;
}

//double f_singlet_integral_test(double T) {
//    double weights_0th[1] = {2.0}; double * w0 = weights_0th;
//    double weights_1st[2] = {1.0, 1.0}; double * w1 = weights_1st;
//    double weights_2nd[3] = {5.0/9.0, 8.0/9.0, 5.0/9.0}; double * w2 = weights_2nd;
//    double weights_3rd[4] = {0.347855, 0.652145, 0.652145, 0.347855}; double * w3 = weights_3rd;
//    double weights_4th[5] = {0.236927, 0.478629, 0.568889, 0.478629, 0.236927}; double * w4 = weights_4th;
//    double *weights[5] = {w0, w1, w2, w3, w4};
//
//    double points_0th[1] = {0}; double *p0 = points_0th;
//    double points_1st[2] = {-1/pow(3,0.5), 1/pow(3,0.5)}; double *p1 = points_1st;
//    double points_2nd[3] = {-pow(3/5,0.5), 0, pow(3/5,0.5)}; double *p2 = points_2nd;
//    double points_3rd[4] = {-0.861136, -0.339981, 0.339981, 0.861136}; double *p3 = points_3rd;
//    double points_4th[5] = {-0.90618, -0.538469, 0, 0.538469, 0.90618}; double *p4 = points_4th;
//
//    double *points[5] = {p0, p1, p2, p3, p4};
//
//    double sum = 0;
//
//    for (int i = 0; i < l; i++) {
//        //double ep = 0 + w_D / (l-1) * i;
//        //sum += 2*f_singlet(ep, T) * w_D / l;
//        sum += w_D * f_singlet(w_D * points[l-1][i], T) * weights[l-1][i];
//    }
//    return sum;
//}
//
//double zero_temp_func(double w, double dE) {
//    if (w == dE)
//        return 0;
//    return 1 / (w - dE);
//}
//
//
//double bound_sign(Vec k, Vec q) {
//    double e_qk = epsilon(k + q) - mu;
//    double e_k = epsilon(k) - mu;
//    int s1 = e_qk < 0;
//    int s2 = e_k < 0;
//    return s1 - s2;
//}
//
//double other_bound_sign(Vec k, Vec q, double T) {
//    double e_qk = epsilon(k + q);
//    double e_k = epsilon(k);
//
//    double f_qk = f(e_qk, T);
//    double f_k = f(e_k, T);
//
//    return f_qk - f_k;
//}
//
//vector<vector<Vec>> get_singularity_freq_surfaces(double (*func)(Vec k, Vec q), Vec q, double center, double width) {
//    vector<vector<Vec>> surfaces;
//    for (int i = 0; i < l; i++) {
//        double s = width * points[l-1][i] + center;
//        vector<Vec> layer = tetrahedron_method(func, q, s);
//        surfaces.push_back(layer);
//    }
//    return surfaces;
//}
//
//vector<vector<Vec>> get_smooth_freq_surfaces(double (*func)(Vec k, Vec q), Vec q, double lower, double upper, double density) {
//    int num_points = int((upper - lower) / density) + 1;
//    vector<vector<Vec>> surfaces;
//    for (int i = 0; i < num_points; i++) {
//        double s = lower + (upper - lower) / (num_points-1) * i;
//        vector<Vec> layer = tetrahedron_method(func, q, s);
//        if(layer.size() > 0)
//            surfaces.push_back(layer);
//    }
//    return surfaces;
//}
//
//double singularity_chi_sum(Vec q, double w, double T, double width) {
//    double de = 0.01;
//
//    double sum = 0;
//    // Get contours of the surface
//    vector<double> surfaces;
//    double upper = w, lower = w - width;
//    for (int both = 0; both < 2; both++) {
//        if (both == 1) {
//            upper = w + width;
//            lower = w;
//        }
//        int num_points = int((upper - lower) / de) + 1;
//        for (int i = 0; i < num_points; i++) {
//            if (both == 1 and i == 0) continue;
//            surfaces.push_back( lower + (upper - lower) / (num_points-1) * i );
//        }
//    }
//
//    for (auto s : surfaces) {
//        vector<Vec> layer = tetrahedron_method(e_diff, q, s);
//        if (s == 0 and w == 0) {
//            sum += get_DOS(layer) * pow(2*M_PI, dim);
//            continue;
//        }
//        for (int j = 0; j < layer.size(); j++) {
//            double r = nonzero_ratio(w, layer[j].freq, layer[j], q, T);
//            if (s == 0 and w == 0) r = 1;
//            sum += r * layer[j].area / vp_diff(layer[j], q);
//        }
//    }
//    return sum * de;
//}
//
//double gauss_chi_sum(Vec q, double w, double T, double width) {
//    double sum = 0;
//    for (int i = 0; i < l; i++) {
//        double s = width * points[l-1][i] + w;
//        vector<Vec> layer = tetrahedron_method(e_diff, q, s);
//        printf("Layer size at s = %f: %d\n", s, layer.size());
//        if (s == 0 and w == 0) {
//            vector<Vec> FS = tetrahedron_method(e_base_avg, q, mu);
//            for (auto k : FS) {
//                if (abs(epsilon(k) - epsilon(k+q)) < 0.01) {
//                    sum += k.area / vp(k) * weights[l-1][i];
//                }
//            }
//            continue;
//        }
//        for (int j = 0; j < layer.size(); j++) {
//            double s_val = width * points[l-1][i] + w;
//            double r = nonzero_ratio(w, s_val, layer[j], q, T) * weights[l-1][i];
//            //if (s == 0 and w == 0) r = weights[l-1][i];
//            sum += r * layer[j].area / vp_diff(layer[j], q);
//            if (r < 0) cout << r << endl;
//        }
//        //printf("Layer %d: %f\n", i, sum);
//    }
//    return sum;
//}
//
//double bound_chi_sum(Vec q, double w, double T, double de, double width, double b, double a) {
//    double sum = 0;
//    // Get contours of the surface
//    vector<double> surfaces;
//    double upper = b, lower = w+width;
//    for (int both = 0; both < 2; both++) {
//        if (both == 1) {
//            upper = w-width;
//            lower = a;
//        }
//        int num_points = int((upper - lower) / de) + 1;
//        for (int i = 0; i < num_points; i++)
//            surfaces.push_back( lower + (upper - lower) / (num_points-1) * i );
//    }
//
//    for (auto s : surfaces) {
//        vector<Vec> layer = tetrahedron_method(e_diff, q, s);
//        for (int j = 0; j < layer.size(); j++) {
//            double r = nonzero_ratio(w, layer[j].freq, layer[j], q, T);
//            sum += r * layer[j].area / vp_diff(layer[j], q);
//        }
//    }
//    return sum * de;
//}
//
//double chi_integrate_freq(Vec q, double width, double w, double T) {
//    double de = 0.01;
//    double a, b;
//    get_bounds(q, w, b, a);
//
//    double sum = 0;
//    sum += gauss_chi_sum(q, w, T, width);
//    //printf("Gauss Sum: %f\n", sum);
//    //sum = singularity_chi_sum(q, w, T, width);
//    //printf("Chi Freq Sum: %f\n", sum);
//    sum += bound_chi_sum(q, w, T, de, width, b, a);
//    //printf("Bound Sum: %f\n", sum);
//    return sum / pow(2*k_max, dim);
//}
//
//void get_bounds(Vec q, double w, double &upper, double &lower) { 
//    vector<Vec> surface = tetrahedron_method(vp_diff, q, 0);
//    double max = 0; double min = 1000;
//    for (int i = 0; i < surface.size(); i++) {
//        double sign = e_diff(surface[i], q) - w;
//        if (sign > 0 and sign > max) max = sign;
//        if (sign < 0 and sign < min) min = sign;
//    }
//    upper = max;
//    lower = min;
//}
//
//double imaginary_chi_integrate(Vec q, double w) {
//    double T = 0.03;
//    double sum = 0;
//    double a, b;
//    vector<Vec> surface1 = tetrahedron_method(e_diff, q, w);
//    //vector<Vec> surface2 = tetrahedron_method(e_diff, q, -w);
//    //cout << "Surface 1: " << surface1.size() << endl;
//    //cout << "Surface 2: " << surface2.size() << endl;
//
//    int counter = 0;
//    //#pragma omp parallel for reduction(+:sum)
//    for (int i = 0; i < surface1.size(); i++) {
//        double sign = other_bound_sign(surface1[i], q, T);
//        sum += sign * surface1[i].area / vp_diff(surface1[i], q);
//        if (sign != 0) counter++;
//        //cout << sign << " " << surface1[i].area << " " << vp_diff(surface1[i], q) << endl;
//        //cout << sum << endl;
//    }
//    //cout << "Sum: " << sum << endl;
//    //cout << "Counter: " << counter << endl;
//
//    //#pragma omp parallel for reduction(+:sum)
//    //for (int i = 0; i < surface2.size(); i++) {
//    //    double sign = other_bound_sign(surface2[i], q);
//    //    sum -= sign * surface2[i].area / vp_diff(surface2[i], q);
//    //    if (sign != 0) counter++;
//    //}
//    //cout << "Counter: " << counter << endl;
//
//    return sum;
//}
//
//void shift_layers(vector<Vec> &layers, Vec shift, double &max, double &min) {
//    max = 0; min = 1000;
//    for (int i = 0; i < layers.size(); i++) {
//        layers[i] = layers[i] + shift;
//        if (epsilon(layers[i]) - mu > max) max = epsilon(layers[i]) - mu;
//        if (epsilon(layers[i]) - mu < min) min = epsilon(layers[i]) - mu;
//    }
//}
//
//double modified_e_diff(Vec k, Vec q) {
//    double eqk = epsilon(k + q), ek = epsilon(k);
//    double f_eqk = f(eqk, 0), f_ek = f(ek, 0);
//    //if (f_eqk == f_ek) return 0;
//    return eqk - ek;
//}
//
//void get_bounds2(Vec q, double &upper, double &lower) {
//    upper = 0; lower = 1000;
//    int pts = 100;
//    for (double i = 0; i < pts; i++) {
//        double x = get_k(i, pts);
//        for (double j = 0; j < pts; j++) {
//            double y = get_k(j, pts);
//            for (double k = 0; k < pts * (dim%2) + 1 * ((dim+1)%2); k++) {
//                double z = get_k(k, pts);
//                Vec k_val(x, y, z);
//                //double f_qk = f(epsilon(k_val + q), 0);
//                //double f_k = f(epsilon(k_val), 0);
//                //if (f_qk != f_k) {
//                    //double diff = epsilon(k_val + q) - epsilon(k_val);
//                    double diff = sphere_func(k_val, q) - sphere_func(k_val, Vec(0,0,0));
//                    if (diff > upper) upper = diff;
//                    if (diff < lower) lower = diff;
//                //}
//            }
//        }
//    }
//}
//
//void get_bounds3(Vec q, double &upper, double &lower) {
//    upper = 0; lower = 1000;
//    int pts = 100;
//    for (double i = 0; i < pts; i++) {
//        double x = get_k(i, pts);
//        for (double j = 0; j < pts; j++) {
//            double y = get_k(j, pts);
//            for (double k = 0; k < pts * (dim%2) + 1 * ((dim+1)%2); k++) {
//                double z = get_k(k, pts);
//                Vec k_val(x, y, z);
//                double val = e_base_avg(k_val, q);
//                //cout << val << endl;
//                if (val > upper) upper = val;
//                if (val < lower) lower = val;
//            }
//        }
//    }
//}
//
//double gauss_chi_sum2(Vec q, double w, double T, double width) {
//    double sum = 0;
//    for (int i = 0; i < l; i++) {
//        double s = width * points[l-1][i] + w;
//        vector<Vec> layer = tetrahedron_method(e_diff, q, s);
//        for (int j = 0; j < layer.size(); j++) {
//            double s_val = width * points[l-1][i] + w;
//            double r = nonzero_ratio(w, s_val, layer[j], q, T) * weights[l-1][i];
//            //if (s == 0 and w == 0) r = weights[l-1][i];
//            sum += r * layer[j].area / vp_diff(layer[j], q);
//            //if (i == 2) cout << r << " " << layer[j].area << " " << vp_diff(layer[j], q) << endl;
//            if (r < 0) cout << r << endl;
//        }
//        printf("Layer %d (size: %d): %f\n", i, layer.size(), sum);
//    }
//    return sum;
//}

//double bound_chi_sum2(Vec q, double w, double T, int pts, double b, double a) {
//    double sum = 0;
//    double A = -4*w + 2*(b+a), B = 4*w - b - 3*a, C = a;
//
//    double prev_s = 0;
//    for (double i = 0; i <= pts; i++) {
//        //double s = A * pow(i/pts, 2) + B * i/pts + C;
//        double s = a + (b-a) * i/pts;
//        if (i == 0) continue;
//        vector<Vec> layer = tetrahedron_method(sphere_func, q, s);
//        for (auto k : layer) {
//            double r = nonzero_ratio(w, k.freq, k, q, T);
//            r = 1;
//            sum += (r * k.area / sphere_func(k, q)) * (s - prev_s);
//        }
//        //cout << "Layer s=" << s << ", ds/dr=" << sphere_func_diff(layer[0],q) 
//        //    << " ds=" << s - prev_s << " area=" << surface_area
//        //    << " (size: " << layer.size() << "): " << sum2 << endl;
//        prev_s = s;
//    }
//    //printf("Sum: %f\n", sum);
//    return sum;
//}
//
//double bound_chi_sum3(Vec q, double w, double T, int pts, double b, double a) {
//    double sum = 0;
//    double prev_s = 0;
//    for (double i = 0; i <= pts; i++) {
//        double s = a + (b-a) * i/pts;
//        if (i == 0) continue;
//        vector<Vec> layer = tetrahedron_method(e_base_avg, q, s);
//        for (auto k : layer) {
//            double r = nonzero_ratio(w, k.freq, k, q, T);
//            sum += r * k.area / vp(k) * (s - prev_s);
//        }
//        prev_s = s;
//    }
//    return sum;
//}

double nonzero_ratio(double w, double dE, Vec k, Vec q, double T) {
    double e_qk = epsilon(k + q) - mu;
    double e_k = epsilon(k) - mu;
    double f_k = f(e_k, T);
    double f_qk = f(e_qk, T);
    if ( fabs(dE - (e_qk - e_k)) > 0.01) {
        cout << "dE: " << dE << " " << e_qk - e_k << endl;
    }

    if (fabs(dE) < 0.0001 and fabs(w) < 0.0001) {
        //return 0;
        if (T == 0 or exp(e_k/T) > 1e6) {
            return e_k < 0;
        }
        double term1 = 1/T * exp(e_k/T) / pow( exp(e_k/T) + 1,2);
        //cout << "Ratio: " << term1 << endl;
        return term1;
    }
    return (f_qk - f_k) / (w - dE);
}

void get_bounds3(Vec q, double &upper, double &lower, double (*func)(Vec k, Vec q)) {
    upper = 0; lower = 1000;
    int pts = 100;
    for (double i = 0; i < pts; i++) {
        double x = get_k(i, pts);
        for (double j = 0; j < pts; j++) {
            double y = get_k(j, pts);
            for (double k = 0; k < pts * (dim%2) + 1 * ((dim+1)%2); k++) {
                double z = get_k(k, pts);
                Vec k_val(x, y, z);
                double val = func(k_val, q);
                if (val > upper) upper = val;
                if (val < lower) lower = val;
            }
        }
    }
}

double sphere_func(Vec k, Vec q) {
    double x = k.vals[0] + q.vals[0], y = k.vals[1] + q.vals[1], z = k.vals[2] + q.vals[2];
    return pow(x,2) + pow(y,2) + pow(z,2);
    //if (k.vals.norm() > 1) return 0;
    return k.vals.squaredNorm();
}

double sphere_func_diff(Vec k, Vec q) {
    double x = k.vals[0] + q.vals[0], y = k.vals[1] + q.vals[1], z = k.vals[2] + q.vals[2];
    return pow(pow(2*x,2) + pow(2*y,2) + pow(2*z,2), 0.5);
    return pow(pow(2*k.vals[0],2) + pow(2*k.vals[1],2) + pow(2*k.vals[2],2),0.5);
}

double denominator(Vec k, Vec q) {
    return k.vals.squaredNorm();
    return (k+q).vals.squaredNorm() - k.vals.squaredNorm();
}

double denominator_diff(Vec k, Vec q) {
    return 2*k.vals.norm();
    Vec temp = 2*(k+q) - 2*k;
    return temp.vals.norm();
}

double integrand_surface(Vec k, Vec q) {
    double x = k.vals[0], y = k.vals[1], z = k.vals[2];
    double r = pow(pow(x,2) + pow(y,2) + pow(z,2), 0.5);
    return r;
    //return (r - 1);
}

double integrand_surface_diff(Vec k, Vec q) {
    double x = k.vals[0], y = k.vals[1], z = k.vals[2];
    double r = pow(pow(x,2) + pow(y,2) + pow(z,2), 0.5);
    return pow(pow(x/r,2) + pow(y/r,2) + pow(z/r,2), 0.5);
}

double integrand(Vec k, Vec q, double w, double T, double (*func)(Vec k, Vec q)) {
    //return 1;
    if (fabs(k.freq - w) < 0.001) return 0;
    return 1 / (k.freq - w);

    //double e_qk = func(k, q) - mu;
    //double e_k = func(k, Vec(0,0,0)) - mu;
    //double f_k = f(e_k, T);
    //double f_qk = f(e_qk, T);

    //double dE = k.freq;
    //if ( fabs(dE - (e_qk - e_k)) > 0.01) {
    //    cout << "dE: " << dE << " " << e_qk - e_k << endl;
    //}

    //if (fabs(dE) < 0.0001 and fabs(w) < 0.0001) {
    //    if (T == 0 or exp(e_k/T) > 1e6)
    //        return e_k < 0;
    //    return 1/T * exp(e_k/T) / pow( exp(e_k/T) + 1,2);
    //}
    //return (f_qk - f_k) / (w - dE);
}

void get_spacing_curve_consts(double w, double a, double b, double &A, double &upr, double &lwr) {
    A = w - a;
    upr = pow((b - w) / (w - a), 1.0/3.0);
    if (isnan(upr)) upr = - pow(fabs(b-w)/fabs(w-a), 1.0/3.0);
    lwr = -1;

    if (w - a > b - w) {
        A = b - w;
        lwr = pow((w - a) / (b - w), 1.0/3.0);
        if (isnan(lwr)) lwr = - pow(fabs(b-w)/fabs(w-a), 1.0/3.0);
        upr = 1;
    }
}

double comparison_integral(Vec q, double w, double b, double a, int pts, double (*func)(Vec k, Vec q)) {
    printf("Starting Comparison Integral\n");
    double sum = 0;
    double A, upr, lwr;
    get_spacing_curve_consts(w, a, b, A, upr, lwr);
    cout << "A: " << A << " Upr: " << upr << " Lwr: " << lwr << endl;
    auto spacing = [A, lwr, upr, w] (double i, double pts) { 
        double x = lwr + (upr- lwr) * i / pts;
        return A * pow(x,3) + w; 
    };

    for (int i = 1; i <= pts; i++) {
        double t = i;
        double r = pow(spacing(t, pts), 0.5);
        double prev_r = pow(spacing(t-1, pts),0.5);
        if (fabs(r*r - w) < 0.001) continue;
        sum += 4 * M_PI * r*r / (r*r - w) * (r - prev_r);
    }
    return sum;
}

double bound_chi_sum4(Vec q, double w, double T, int pts, double b, double a, double (*func)(Vec k, Vec q), double (*func_diff)(Vec k, Vec q)) {
    double sum = 0;

    double A, upr, lwr;
    get_spacing_curve_consts(w, a, b, A, upr, lwr);
    auto spacing = [A, w, lwr, upr] (double i, double pts) { 
        double x = lwr + (upr- lwr) * i / pts;
        return A * pow(x,3) + w; 
    };

    #pragma omp parallel for reduction(+:sum)
    for (int i = 1; i <= pts; i++) {
        double t = i;
        double s = spacing(t, pts);
        double prev_s = spacing(t-1, pts);

        vector<Vec> layer = tetrahedron_method(func, q, s);
        for (auto k : layer) {
            double r = integrand(k, q, w, T, func);
            sum += r * k.area / func_diff(k,q) * fabs(s - prev_s);
        }
    }
    return sum;
}


double chi_ep_integrate(Vec q, double w, double T) {
    double a, b;
    get_bounds3(q, b, a, denominator);
    b = 4;

    double sum = 0;
    sum += bound_chi_sum4(q, w, T, 100, b, a, denominator, denominator_diff);
    return sum; // / pow(2*k_max, dim);
}
