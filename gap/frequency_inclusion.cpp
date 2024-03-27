#include <iostream>
#include <unistd.h>
#include <ctime>

#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#include <algorithm>
#include <unordered_map>
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
        vector<Vec> layer = tetrahedron_method(ep);
        basis.push_back(layer);
    }
    return basis;
}

MatrixXd create_P_freq(vector<vector<Vec>> &k, double T, const vector<vector<vector<double>>> &chi_cube2) {
    int size = 0;
    for (int i = 0; i < k.size(); i++)
        size += k[i].size();
    double avg = 1.0 * size / k.size();

    MatrixXd P(size, size);

    double weights_0th[1] = {2.0}; double * w0 = weights_0th;
    double weights_1st[2] = {1.0, 1.0}; double * w1 = weights_1st;
    double weights_2nd[3] = {5.0/9.0, 8.0/9.0, 5.0/9.0}; double * w2 = weights_2nd;
    double weights_3rd[4] = {0.347855, 0.652145, 0.652145, 0.347855}; double * w3 = weights_3rd;
    double weights_4th[5] = {0.236927, 0.478629, 0.568889, 0.478629, 0.236927}; double * w4 = weights_4th;
    double *weights[5] = {w0, w1, w2, w3, w4};

    double points_0th[1] = {0}; double *p0 = points_0th;
    double points_1st[2] = {-1/pow(3,0.5), 1/pow(3,0.5)}; double *p1 = points_1st;
    double points_2nd[3] = {-pow(3/5,0.5), 0, pow(3/5,0.5)}; double *p2 = points_2nd;
    double points_3rd[4] = {-0.861136, -0.339981, 0.339981, 0.861136}; double *p3 = points_3rd;
    double points_4th[5] = {-0.90618, -0.538469, 0, 0.538469, 0.90618}; double *p4 = points_4th;

    double *points[5] = {p0, p1, p2, p3, p4};

    cout << "Creating P Matrix with frequency\n";
    int ind1 = 0;
    for (int i = 0; i < k.size(); i++) {
        #pragma omp parallel for
        for (int j = 0; j < k[i].size(); j++) {
            Vec k1 = k[i][j];
            int ind2 = 0;
            for (int x = 0; x < k.size(); x++) {
                for (int y = 0; y < k[x].size(); y++) {
                    Vec k2 = k[x][y];
                    double d1 = pow(k1.area/vp(k1),0.5); 
                    double d2 = pow(k2.area/vp(k2),0.5); 
                    double f1 = f_singlet(w_D * points[l-1][x], T);
                    // f * d_epsilon
                    double fde = f_singlet(w_D * points[l-1][x], T) * weights[l-1][x];

                    P(ind1,ind2) = - d1 * d2 * fde * V(k1, k2, T, chi_cube2); 
                    //P(ind1, ind2) = f_singlet(w_D*points[l-1][x], T) * weights[l-1][x];
                    ind2++;
                }
            }
            ind1++;
        }
    }
    cout << "P Matrix Created\n";

    //return P * w_D / (k[0].size());
    return P * w_D * (2 / pow(2*M_PI, dim)); 
}

MatrixXd create_P_freq2(vector<vector<Vec>> &k, double T, const vector<vector<vector<double>>> &chi_cube2) {
    int size = 0;
    for (int i = 0; i < k.size(); i++) {
        size += k[i].size();
    }

    MatrixXd P(size, size);
    cout << "Creating P Matrix with frequency\n";
    int ind1 = 0;
    for (int i = 0; i < k.size(); i++) {
        for (int j = 0; j < k[i].size(); j++) {
            Vec k1 = k[i][j];
            #pragma omp parallel for
            for (int x = 0; x < k.size(); x++) {

                int ind2 = 0;
                for (int temp = 0; temp < x; temp++)
                    ind2 += k[temp].size();

                for (int y = 0; y < k[x].size(); y++) {
                    Vec k2 = k[x][y];
                    double d1 = pow(k1.area/vp(k1),0.5); 
                    double d2 = pow(k2.area/vp(k2),0.5); 
                    double f1 = f_singlet(w_D * points[l-1][x], T);
                    // f * d_epsilon
                    double fde1 = f_singlet(w_D * points[l-1][i], T) * weights[l-1][i];
                    double fde2 = f_singlet(w_D * points[l-1][x], T) * weights[l-1][x];

                    P(ind1,ind2) = - d1 * d2 * pow(fde1*fde2,0.5) * V(k1, k2, T, chi_cube2); 
                    ind2++;
                }
            }
            ind1++;
        }
    }
    cout << "P Matrix Created\n";

    //return P * 2 * w_D / (l * k_size);
    return P * w_D * (2 / pow(2*M_PI, dim)); 
}

double f_singlet_integral_test(double T) {
    double weights_0th[1] = {2.0}; double * w0 = weights_0th;
    double weights_1st[2] = {1.0, 1.0}; double * w1 = weights_1st;
    double weights_2nd[3] = {5.0/9.0, 8.0/9.0, 5.0/9.0}; double * w2 = weights_2nd;
    double weights_3rd[4] = {0.347855, 0.652145, 0.652145, 0.347855}; double * w3 = weights_3rd;
    double weights_4th[5] = {0.236927, 0.478629, 0.568889, 0.478629, 0.236927}; double * w4 = weights_4th;
    double *weights[5] = {w0, w1, w2, w3, w4};

    double points_0th[1] = {0}; double *p0 = points_0th;
    double points_1st[2] = {-1/pow(3,0.5), 1/pow(3,0.5)}; double *p1 = points_1st;
    double points_2nd[3] = {-pow(3/5,0.5), 0, pow(3/5,0.5)}; double *p2 = points_2nd;
    double points_3rd[4] = {-0.861136, -0.339981, 0.339981, 0.861136}; double *p3 = points_3rd;
    double points_4th[5] = {-0.90618, -0.538469, 0, 0.538469, 0.90618}; double *p4 = points_4th;

    double *points[5] = {p0, p1, p2, p3, p4};

    double sum = 0;

    for (int i = 0; i < l; i++) {
        //double ep = 0 + w_D / (l-1) * i;
        //sum += 2*f_singlet(ep, T) * w_D / l;
        sum += w_D * f_singlet(w_D * points[l-1][i], T) * weights[l-1][i];
    }
    return sum;
}

