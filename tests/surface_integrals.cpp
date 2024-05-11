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

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::unordered_map;
//using namespace Eigen;

double trapezoid_area_from_points(Vec k1, Vec k2, Vec k3, Vec k4) {
    // Define distances 
    Vec k12 = k1 - k2; if (k12.cartesian) k12.to_spherical();
    Vec k23 = k2 - k3; if (k23.cartesian) k23.to_spherical();
    Vec k34 = k3 - k4; if (k34.cartesian) k34.to_spherical();
    Vec k14 = k1 - k4; if (k14.cartesian) k14.to_spherical();
    Vec k13 = k1 - k3; if (k13.cartesian) k13.to_spherical();
    double d12 = k12.vals(0);
    double d23 = k23.vals(0);
    double d34 = k34.vals(0);
    double d14 = k14.vals(0);
    double d13 = k13.vals(0);

    double s = (d12 + d23 + d34 + d14)/2; 
    double theta1 = acos( (pow(d12,2) + pow(d23,2) - pow(d13,2)) / (2*d12*d23));
    double theta2 = acos( (pow(d14,2) + pow(d34,2) - pow(d13,2)) / (2*d14*d34));
    double theta = theta1 + theta2;
    double area = pow((s-d12)*(s-d23)*(s-d34)*(s-d14) - d12*d23*d34*d14*pow(cos(theta/2),2), 0.5);

    return area;
}

vector<double> get_wave_from_file() {
    ifstream file("../data/const3D_mu=-1.2_U=4.0_wD=1.0_n=6.dat");
    string line;
    getline(file, line);
    vector<double> wave;
    while (file) {
        double x, y, z, c1;
        file >> x >> y >> z >> c1;
        wave.push_back(c1);
        string line;
        getline(file, line);
    }
    return wave;
}

double test_surface_integral_discrete(vector<Vec> &FS_vecs, int i, double T, const vector<vector<vector<double>>> &chi_cube, auto &V, bool norm = false) {
    auto ind_func = [&V](auto &FS_vecs, auto &dumb_wave, int i, int j, double T, const vector<vector<vector<double>>> &chi_cube, bool norm = false) {
        Vec k1 = FS_vecs[i]; Vec k2 = FS_vecs[j];
        if (not norm) {
            return dumb_wave[i] * V(k1, k2, T, chi_cube) * dumb_wave[j] / (vp(k1) * vp(k2));
        }
        return dumb_wave[j] * dumb_wave[j];// / (vp(k2)*vp(k2));
    };
    vector<double> dumb_wave = get_wave_from_file();
    double sum = 0;
    for (int j = 0; j < FS_vecs.size(); j++) {
        Vec k1 = FS_vecs[j];
        double f = ind_func(FS_vecs, dumb_wave, i, j, T, chi_cube, norm);
        sum += f*k1.area;
    }
    return sum;
}

double test_4d_surface_integral_discrete(vector<Vec> &FS_vecs, const vector<vector<vector<double>>> &chi_cube, double T, auto &V) {
    
    double normalized = test_surface_integral_discrete(FS_vecs, 0, T, chi_cube, V, true);
    cout << "Discrete Normalization: " << normalized << endl;
    double sum = 0;
    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < FS_vecs.size(); i++) {
        Vec k = FS_vecs[i];
        double f = test_surface_integral_discrete(FS_vecs, i, T, chi_cube, V);
        sum += f*k.area;
    }
    return -sum / normalized;

}

double test_surface_integral4(vector<Vec> &FS_vecs, Vec p, double T, double mu, const vector<vector<vector<double>>> &chi_cube, auto &func) {
    double sum = 0;
    for (int i = 0; i < FS_vecs.size(); i++) {
        Vec k1 = FS_vecs[i];
        double f = func(k1, p, T, mu, chi_cube);
        sum += f*k1.area;
    }
    return sum ;//* (2 / pow(2*M_PI,3));
}

double test_4d_surface_integral4(vector<Vec> &FS_vecs, const vector<vector<vector<double>>> &chi_cube, double T, double mu, auto &wave, auto &V) {

    auto norm_func = [&wave](Vec k, Vec extra, double T, double mu, const vector<vector<vector<double>>> &chi_cube) {
        return pow(wave(k, mu),2) / vp(k);
    };
    auto func = [&wave, &V](Vec p1, Vec p2, double T, double mu, const vector<vector<vector<double>>> &chi_cube) { 
        //return V(p1, p2, T, chi_cube);
        return wave(p1, mu) * V(p1,p2,T,chi_cube) * wave(p2, mu) / (vp(p1) * vp(p2)); 
    };

    double normalized = test_surface_integral4(FS_vecs, FS_vecs[0], T, mu, chi_cube, norm_func);
    //cout << "Normalized: " << normalized << endl;

    double sum = 0;
    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < FS_vecs.size(); i++) {
        Vec k = FS_vecs[i];
        double f = test_surface_integral4(FS_vecs, k, T, mu, chi_cube, func);
        sum += f*k.area;
    }
    //cout << "Sum: " << sum << endl;
    return -sum / normalized * (2 / pow(2*M_PI,dim));
}


void test_integral_pairing_interaction_parity(vector<Vec> &FS, double T, double mu) {

    auto wave = [](Vec k, double mu) { 
        Vec q = k; 
        if (q.cartesian == false) q.to_cartesian();
        return cos(q.vals(1))-cos(q.vals(0)); 
    };
    auto Vs = [](Vec p1, Vec p2, double T, const auto &chi_cube) { 
        return V(p1, p2, 0, T, chi_cube);
    };
    auto func = [&wave, &Vs](Vec p1, Vec p2, double T, double mu, const vector<vector<vector<double>>> &chi_cube) { 
        return wave(p1, mu) * Vs(p1,p2,T,chi_cube) * wave(p2, mu);// / (vp(p1) * vp(p2)); 
    };

    double DOS = 0; for (auto x : FS) DOS += x.area / vp(x);
    DOS /= pow(2*M_PI,3);
    auto cube = chi_cube(T, mu, DOS, 0);
    double num = 10;
    for (double i = 0; i < num; i++) {
        for (double j = 0; j < num; j++) {
            for (double k = 0; k < num; k++) {
                double x = 2*k_max*(i/num-1.0);
                double y = 2*k_max*(j/num-1.0);
                double z = 2*k_max*(k/num-1.0);

                Vec p(x, y, z);

                double val1 = test_surface_integral4(FS, p, T, mu, cube, func);
                double val2 = test_surface_integral4(FS, -1*p, T, mu, cube, func);
                if ( fabs(val1 - val2) > 0.0000001) {
                    cout << setprecision(10);
                    cout << "P: " << p << endl;
                    cout << "Integral: " << val1 << " " << val2 << endl;
                    cout << "Chi: " << calculate_chi_from_cube(cube, p) << endl;
                }
            }
        }
    }

}

double heavy_step(Vec k) {
    return k.vals.norm() < 1;
}

double ratio_1D(Vec q, Vec k) {
    double e_k = epsilon_sphere(k);
    double e_qk = epsilon_sphere(k+q);

    double f_k = heavy_step(k);
    double f_qk = heavy_step(k+q);

    return (f_qk - f_k) / (e_k - e_qk);
}

double explicit_formula_1d(Vec q, Vec k) {
    return 1 / (q.vals.squaredNorm() - 4*k.vals.squaredNorm());
}

double example_formula_1D(Vec q, Vec k) {
    //return (k - q).vals.norm();
    return 1/(k.vals(0) - q.vals(0));
}

bool k_in_dense_box(Vec k, vector<Vec> singularities, vector<Vec> singularity_widths) {
    for (int i = 0; i < singularities.size(); i++) {
        Vec singularity = singularities[i];
        Vec singularity_width = singularity_widths[i];
        if (k.vals(0) >= singularity.vals(0) - singularity_width.vals(0) 
            and k.vals(0) <= singularity.vals(0) + singularity_width.vals(0)
            and k.vals(1) >= singularity.vals(1) - singularity_width.vals(1) 
            and k.vals(1) <= singularity.vals(1) + singularity_width.vals(1) 
            and k.vals(2) >= singularity.vals(2) - singularity_width.vals(2) 
            and k.vals(2) <= singularity.vals(2) + singularity_width.vals(2)) {
            return true;
        }
    }
    return false;
}

void plot_coupling() {

    auto dx2_y2 = [](Vec k, double mu) { 
        Vec q = k; 
        if (q.cartesian == false) q.to_cartesian();
        return cos(q.vals(1))-cos(q.vals(0)); 
    };
    auto dz2_1 = [](Vec k, double mu) { 
        Vec q = k; 
        if (q.cartesian == false) q.to_cartesian();
        return 2*cos(q.vals(2))-cos(q.vals(0))-cos(q.vals(1)); 
    };
    auto dxy = [](Vec k, double mu) { 
        Vec q = k; if (not q.cartesian) q.to_cartesian();
        return sin(q.vals(0))*sin(q.vals(1));
    };
    auto dxz = [](Vec k, double mu) { 
        Vec q = k; if (not q.cartesian) q.to_cartesian();
        return sin(q.vals(0))*sin(q.vals(2));
    };
    auto dyz = [](Vec k, double mu) { 
        Vec q = k; if (not q.cartesian) q.to_cartesian();
        return sin(q.vals(1))*sin(q.vals(2));
    };
    auto px = [](Vec k, double mu) { 
        Vec q = k; if (not q.cartesian) q.to_cartesian();
        return sin(q.vals(0));
    };
    auto py = [](Vec k, double mu) { 
        Vec q = k; if (not q.cartesian) q.to_cartesian();
        return sin(q.vals(1));
    };
    auto pz = [](Vec k, double mu) { 
        Vec q = k; if (not q.cartesian) q.to_cartesian();
        return sin(q.vals(2));
    };
    auto p_wave2 = [](Vec k, double mu) { 
        Vec q = k; if (not q.cartesian) q.to_cartesian();
        return (sin(q.vals(0)) - sin(q.vals(1)));
    };
    auto extended_s_wave = [](Vec k, double mu) { 
        return -mu/2 + cos(k.vals(0)) + cos(k.vals(1)) + cos(k.vals(2)); 
    };
    auto wave_norm = [](Vec k, double mu) {
        return 1.0; 
    };
    auto Vs = [](Vec p1, Vec p2, double T, const auto &chi_cube) { 
        return V(p1, p2, T, chi_cube);
    };
    auto Vt = [](Vec p1, Vec p2, double T, const auto &chi_cube) { 
        Vec q = to_IBZ_2(p1 - p2);
        double chi = calculate_chi_from_cube(chi_cube, q);
        return -pow(U,2)*chi / (1 - pow(U*chi,2));
    };
    auto V_norm = [](Vec p1, Vec p2, double T, const auto &chi_cube) { 
        Vec q = to_IBZ_2(p1 - p2);
        double chi = calculate_chi_from_cube(chi_cube, q);
        double V = pow(U,3)*pow(chi,2) / (1-U*chi) + pow(U,2)*chi / (1-pow(U*chi,2));
        return -V;
    };

    double T = 0.0065;
    T = 0.25;
    double mu = -0.9;
    ofstream file("coupling.dat");

    while (mu > -4.5) {
        cout << "Mu: " << mu << endl;
        vector<Vec> FS = tetrahedron_method(mu); 
        cout << "Fermi Surface Size: " << FS.size() << endl;
        double DOS = get_DOS(FS);

        vector<vector<vector<double>>> cube = chi_cube(T, mu, DOS, 0);
        MatrixXd P = create_P(FS, T, cube);
        vector<EigAndVec> solutions = power_iteration(P, 0.0001);
        cout << "Power Iteration Done\n";
        //cout << P.eigenvalues().real().maxCoeff();

        double eig = solutions[solutions.size()-1].eig;
        //double eig = 0;
        cout << "Eig Found: " << eig << " \n";

        // Renormalization
        double renormalization = test_4d_surface_integral4(FS, cube, T, mu, wave_norm, V_norm)+1.0;
        //double renormalization = 1;
        //cout << "\nEigs: " << endl;
        //for (auto x : solutions) {
        //    cout << x.eig / renormalization << endl;
        //}
        //cout << endl;
        // D-waves
        double dx2_y2_coupling = test_4d_surface_integral4(FS, cube, T, mu, dx2_y2, Vs);
        //double dz2_1_coupling, dxy_coupling, dxz_coupling, dyz_coupling, px_coupling, py_coupling, pz_coupling, s_coupling, extended_s_coupling, dx_coupling = 0;
        //double dz2_1_coupling = test_4d_surface_integral4(FS, cube, T, mu, dz2_1, Vs);
        double dxy_coupling = test_4d_surface_integral4(FS, cube, T, mu, dxy, Vs);
        //double dxz_coupling = test_4d_surface_integral4(FS, cube, T, mu, dxz, Vs);
        //double dyz_coupling = test_4d_surface_integral4(FS, cube, T, mu, dyz, Vs);
        // P-waves
        double px_coupling = test_4d_surface_integral4(FS, cube, T, mu, px, Vt);
        //double py_coupling = test_4d_surface_integral4(FS, cube, T, mu, py, Vt);
        //double pz_coupling = test_4d_surface_integral4(FS, cube, T, mu, pz, Vt);

        //// P-wave for singlet potential
        //double dx_coupling = test_4d_surface_integral4(FS, cube, T, mu, px, Vs);
        ////double dy_coupling = test_4d_surface_integral4(FS, cube, T, mu, py, Vs);
        ////double dz_coupling = test_4d_surface_integral4(FS, cube, T, mu, pz, Vs);
        ////double p_coupling2 = test_4d_surface_integral4(FS, cube, T, mu, p_wave2, Vs);
        //double s_coupling = test_4d_surface_integral4(FS, cube, T, mu, wave_norm, Vs);
        double extended_s_coupling = test_4d_surface_integral4(FS, cube, T, mu, extended_s_wave, Vs);
        //double dumb_wave_coupling = test_4d_surface_integral_discrete(FS, cube, T, Vs);
        //double dz2_1_coupling = 0, dxz_coupling = 0, dyz_coupling = 0, px_coupling = 0, py_coupling = 0, pz_coupling = 0, s_coupling = 0; 
        file << setprecision(4);
        file << mu << " ";
        file << renormalization << " ";

        file << dx2_y2_coupling / renormalization << " "<< " ";
        //file << dz2_1_coupling / renormalization << " "<< " ";
        file << dxy_coupling / renormalization << " "<< " ";
        //file << dxz_coupling / renormalization << " "<< " ";
        //file << dyz_coupling / renormalization << " "<< " ";

        file << px_coupling / renormalization << " "<< " ";
        //file << py_coupling / renormalization << " "<< " ";
        //file << pz_coupling / renormalization << " "<< " ";

        //cout << dy_coupling / renormalization << " ";
        //cout << dz_coupling / renormalization << " ";
        //cout << "P-coupling2:  p_coupling2 / renormalization << " ";
        //cout << s_coupling << " " << s_coupling / renormalization << endl;
        //file << s_coupling / renormalization << " ";
        file << extended_s_coupling / renormalization  << " ";
        //cout << "Dumb Wave:  dumb_wave_coupling / renormalization  << " ";

        //double eig = 0;
        file << eig / renormalization << " " << " ";
        file << endl;
        mu -= 0.1;
    }

}

void plot_DOS() {
    ofstream file("info.log");
    double num = 100;
    vector<double> DOS(num);
    double mu_max = 5.9;
    for (int i = 0; i < num; i++) {
        double new_mu = -mu_max + 2.0*mu_max*i/(num-1);
        new_mu = 0.0;
        init_config(mu, U, t, tn, w_D, new_mu, U, t, tn, w_D);

        auto wave = [](Vec k, double mu) { 
            return 1.0;
        };
        auto norm_func = [&wave](Vec k, Vec extra, double T, double mu, const vector<vector<vector<double>>> &chi_cube) {
            //return 1.0;
            return pow(wave(k, mu),2) / vp(k);
        };

        vector<Vec> FS = tetrahedron_method(mu);
        filter_FS(FS);
        double T = 0.0065;
        vector<vector<vector<double>>> cube;
        double test = test_surface_integral4(FS, FS[0], T, mu, cube, norm_func);
        DOS[i] = test;
        cout << new_mu << " " << DOS[i] << endl;
        break;
    }
}
