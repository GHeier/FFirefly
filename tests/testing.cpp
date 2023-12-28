#include <iomanip>
#include <iostream>
#include <fstream>
//#include <chrono>

#include <cstring>
#include <string>
#include <math.h>
#include <vector>

#include <unordered_set>
#include <unordered_map>
#include <boost/functional/hash.hpp>
#include <gsl/gsl_integration.h>

//#include <boost/math/tools/roots.hpp>
//#include <Eigen/Dense>
//#include <gsl/gsl_integration.h>
//#include <python3.10/Python.h>

#include "cfg.h"
#include "fermi_surface.h"
#include "vec.h"
#include "analysis.h"
#include "band_structure.h"
//#include "py_port.h"
#include "save_data.h"
#include "utilities.h"
#include "potential.h"

#include "plot.cpp"
#include "surface_integrals.cpp"


using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::unordered_map;
//using namespace Eigen;


bool test_coordinate_transforms() {
    Vec x(1.0, 0.0, 0.0, false);
    Vec y(1.5, 1.0, 0.0, false);
    Vec z(2.0, 0.0, 0.0, false);
    Vec a(1.0, 0.0, M_PI/2, false);
    Vec b(1.0, 1.4, -0.5, false);
    Vec c(6.06, 0.34, -0.96, false);

    Vec x1 = x; Vec y1 = y; Vec z1 = z; Vec a1 = a; Vec b1 = b; Vec c1 = c;

    x.to_cartesian();
    y.to_cartesian();
    z.to_cartesian();
    a.to_cartesian();
    b.to_cartesian();
    c.to_cartesian();

    if (fabs(x.vals(0) - 1.0) > 0.00001 or 
            fabs(x.vals(1) - 0.0) > 0.00001)
        cout << "Transform 1 Failed\n";
    if (fabs(y.vals(0) - 0.810453) > 0.00001 or 
            fabs(y.vals(1) - 1.26221) > 0.00001) 
        cout << "Transform 2 Failed\n";
    if (fabs(z.vals(0) - 2) > 0.00001 or 
            fabs(z.vals(1) - 0.0) > 0.00001 or 
            fabs(z.vals(2) - 0.0) > 0.00001) 
        cout << "Transform 3 Failed\n";
    if (fabs(a.vals(0) - 0) > 0.00001 or 
            fabs(a.vals(1) - 0.0) > 0.00001 or 
            fabs(a.vals(2) - 1.0) > 0.00001) 
        cout << "Transform 4 Failed\n";
    if (fabs(b.vals(0) - 0.149160) > 0.00001 or 
            fabs(b.vals(1) - 0.864813) > 0.00001 or 
            fabs(b.vals(2) - (-0.479426)) > 0.00001) 
        cout << "Transform 5 Failed\n";
    if (fabs(c.vals(0) - 3.27657) > 0.0001 or 
            fabs(c.vals(1) - 1.15904) > 0.0001 or 
            fabs(c.vals(2) - (-4.9643)) > 0.0001) 
        cout << "Transform 6 Failed\n";

    x.to_spherical();
    y.to_spherical();
    z.to_spherical();
    a.to_spherical();
    b.to_spherical();
    c.to_spherical();

    for (int i = 0; i < dim; i++) if (fabs(x.vals(i) - x1.vals(i)) > 0.00001) 
        cout << "Cartesian to spherical Transform #1 Failed";
    for (int i = 0; i < dim; i++) if (fabs(y.vals(i) - y1.vals(i)) > 0.00001) 
        cout << "Cartesian to spherical Transform #2 Failed";
    for (int i = 0; i < dim; i++) if (fabs(z.vals(i) - z1.vals(i)) > 0.00001) 
        cout << "Cartesian to spherical Transform #3 Failed";
    for (int i = 0; i < dim; i++) if (fabs(a.vals(i) - a1.vals(i)) > 0.00001) 
        cout <<"Cartesian to spherical Transform #4 Failed";
    for (int i = 0; i < dim; i++) if (fabs(b.vals(i) - b1.vals(i)) > 0.00001) 
        cout << "Cartesian to spherical Transform #5 Failed";
    for (int i = 0; i < dim; i++) if (fabs(c.vals(i) - c1.vals(i)) > 0.00001) 
        cout << "Cartesian to spherical Transform #6 Failed";
    cout << "All Transforms Completed. If no error messages, the test is passed\n";

    return true;
}

void test_get_k() {
    int n = 5;
    double correct_values[n] = {-M_PI, -M_PI/2, 0, M_PI/2, M_PI};
    for (int i = 0; i < n; i++) {
        double value = get_k(i, n);
        if ( fabs(value - correct_values[i]) > 0.000001)
            cout << "get_k() failed test.\n";
    }
}

void test_set(vector<Vec> FS_vecs) {
    vector<string> q;
    for (int i = 0; i < FS_vecs.size(); i++) {
        for (int j = 0; j < FS_vecs.size(); j++) {
            Vec k_plus = to_IBZ_spherical(FS_vecs[i] + FS_vecs[j]);
            Vec k_minus = to_IBZ_spherical(FS_vecs[i] - FS_vecs[j]);
            q.push_back(vec_to_string(k_plus));
            q.push_back(vec_to_string(k_minus));
        }
    }
    unordered_set<string> s( q.begin(), q.end() );
    for (auto x: s)
        cout << x << endl;
    cout << q.size() << endl;
    cout << s.size() << endl;

}

void test_IBZ2() {
    int num = 10;
    for (double i = 0; i < num; i++) {
        for (double j = 0; j < num; j++) {
            for (double k = 0; k < num; k++) {
                double x = get_k(i, num);
                double y = get_k(j, num);
                double z = get_k(k, num);
                Vec input(x, y, z);
                Vec output = to_IBZ_2(input);
                double chi1 = chi_trapezoidal(input, 0.25, -1.2, 40);
                double chi2 = chi_trapezoidal(output, 0.25, -1.2, 40);
                if (fabs(chi1 - chi2) > 0.0000001) 
                    cout << "Failed for " << input << "->" << output << endl;
            }
        }
    }
}

void test_IBZ() {
    double T = 0.25, mu = -1.2;
    Vec q1(1, 2, .7);
    Vec q2(.1, -2, .7);
    Vec q3(-.8, -1.2, -1.1);
    double chi1 = chi_trapezoidal(q1, T, mu, 30);
    double chi2 = chi_trapezoidal(q2, T, mu, 30);
    double chi3 = chi_trapezoidal(q3, T, mu, 30);
    Vec q1_t = to_IBZ_2(q1);
    Vec q2_t = to_IBZ_2(q2);
    Vec q3_t = to_IBZ_2(q3);
    double chi1_t = chi_trapezoidal(q1_t, T, mu, 30);
    double chi2_t = chi_trapezoidal(q2_t, T, mu, 30);
    double chi3_t = chi_trapezoidal(q3_t, T, mu, 30);
    cout << chi1 << " " << chi1_t << endl;
    cout << chi2 << " " << chi2_t << endl;
    cout << chi3 << " " << chi3_t << endl;
}

void test_diagonalization() {
    int n = 20;
    auto func = [](double x, double y, double a, double b) {
        return cos(x)*cos(y) + 0.0;
        return (cos(x)-cos(y))*(cos(a)-cos(b)) + 1.0;
    };
    MatrixXd V(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double theta1 = 6.28*i/n, phi1 = 3.14*j/n;
            double x = cos(theta1)*sin(phi1), y = sin(theta1)*sin(phi1);
            V(i,j) = func(theta1, 2*phi1, 0, 0);
            //for (int k = 0; k < n; k++) {
            //    for (int l = 0; l < n; l++) {
            //        double theta1 = 6.28*i/n, phi1 = 3.14*j/n;
            //        double theta2 = 6.28*k/n, phi2 = 3.14*l/n;
            //        double x = cos(theta1)*sin(phi1), y = sin(theta1)*sin(phi1);
            //        double a = cos(theta2)*sin(phi2), b = sin(theta2)*sin(phi2);
            //        int ind1 = i*n + j, ind2 = k*n + l;
            //        V(ind1,ind2) = func(x, y, a, b);
            //        //V(ind1,ind2) = func(theta1, phi1, theta2, phi2);
            //    }
            //}
        }
    }
    EigenSolver<MatrixXd> s(V);
    VectorXcd vals = s.eigenvalues() / (n);
    EigenSolver<MatrixXd>::EigenvectorsType vecs = s.eigenvectors();

    vector<EigAndVec> solutions = combine_eigs_and_vecs(vals.real(), vecs.real());
    sort(solutions.rbegin(), solutions.rend());

    ofstream file("1val.dat");
    ofstream file2("2val.dat");
    cout << solutions[0].eig << endl << solutions[1].eig << endl;
    file << solutions[0].vec;
    file2 << solutions[1].vec;
}

void test_pairing_interaction_parity(vector<Vec> &FS, double T, double mu) {
    double DOS = 0; for (auto x : FS) DOS += x.area / vp(x);
    DOS /= pow(2*M_PI, 3);
    auto cube = chi_cube(T, mu, DOS);
    Vec k2;
    double num = 10;
    for (double i = 0; i < num; i++) {
        for (double j = 0; j < num; j++) {
            for (double k = 0; k < num; k++) {
                double x = get_k(i, num);
                double y = get_k(j, num);
                double z = get_k(k, num);

                Vec k1(x, y, z);

                double V1 = potential_scalapino_cube(k1, k2, T, cube);
                double V2 = potential_scalapino_cube(-1*k1, k2, T, cube);

                if ( fabs(V1 - V2) > 0.00001) 
                    cout << k1 << "->" << V1 << ", " << V2 << endl;
            }
        }
    }

}

int f(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) {
    //double sigma = *((double *) fdata); // we can pass σ via fdata argument
    Vec q = *((Vec *) fdata);
    //double sum = 0;
    //unsigned i;
    //for (i = 0; i < ndim; ++i) sum += x[i] * x[i];
    //// compute the output value: note that fdim should == 1 from below
    //fval[0] = exp(-sigma * sum);
    Vec k_val(x[0], x[1], x[2]);
    double T = 0.25;
    fval[0] = ratio(q, k_val, T, mu);
    //cout << k_val << fval[0] << endl;
    //fval[0] = x[0]*x[0] * x[1]*x[1] * x[2]*x[2];
    return 0; // success*
}


int main() {
    //plot_chi(0.25);
    //plot_chi2(0.25);
    plot_coupling();
    return 0;
}
