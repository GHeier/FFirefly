#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <string>
#include <vector>
#include <algorithm>

#include <omp.h>
#include <boost/functional/hash.hpp>
#include <tuple>
#include <unordered_map>

#include "utilities.h"
#include "potential.h"
#include "vec.h"
#include "fermi_surface.h"
#include "band_structure.h"
#include "cfg.h"
#include "frequency_inclusion.hpp"

using namespace std;


double potential_const(Vec k1, Vec k2) {
    return -1;
}

double potential_test(Vec k1, Vec k2) {
    double volume_scaling = (2 / pow(2*M_PI, dim));
    return -1.0;
    Vec q1 = k1;
    Vec q2 = k2;
    if (q1.cartesian == false) q1.to_cartesian();
    if (q2.cartesian == false) q2.to_cartesian();
    //return -1 * pow(sin(q1.vals[0])*sin(q2.vals[0]),2);
    //return -1 * (sin(q1.vals[0]) - sin(q1.vals[1]))*(sin(q2.vals[0]) - sin(q2.vals[1]));
    //return -1*(sin(q1.vals[0])*sin(q1.vals[1]) * sin(q2.vals[0])*sin(q2.vals[1]));
    //return ( cos( q1.vals[0] + q1.vals[1] ) )*( cos( q2.vals[0] + q2.vals[1] ) );
    //return ( cos( q1.vals[0] ) - cos( q1.vals[1] ) )*( cos( q2.vals[0] ) + cos( q2.vals[1] ) );
    //return cos(q1.vals[0])*cos(q2.vals[0])/M_PI + 1 / (2*M_PI);
    //return -1*( sin(q1.vals[0]) + sin(q1.vals[1]) )*( sin(q2.vals[0]) + sin(q2.vals[1]) );
    return -1*( cos(q1.vals[0]) - cos(q1.vals[1]) )*( cos(q2.vals[0]) - cos(q2.vals[1]) ) + (-0.5)*sin(q1.vals[0])*sin(q1.vals[1])*sin(q2.vals[0])*sin(q2.vals[1]);
    //return 100*( 2*cos(q1.vals[2]) - cos( q1.vals[0] ) - cos( q1.vals[1] ) ) 
    //    * ( 2*cos(q2.vals[2]) - cos( q2.vals[0] ) - cos( q2.vals[1] ) ) + 0;
}
double phonon_coulomb(Vec q) {
    if (q.cartesian == false) q.to_cartesian();
    double qx = q.vals[0];
    double Vp = 1/3;
    if (q.norm() != 0) {
        Vp = 1/(1+2*qx*qx / pow(q.norm(), 2));
    }
    double Vc = 1 / (1 + q.norm());
    return Vp + Vc;
}

double potential_scal(Vec k1, Vec k2, double T) {
    Vec q2 = k1 - k2;
    Vec q = to_IBZ_2(q2);
    
    double chi_sub = chi_trapezoidal(q, T, mu, 0, 30);
    if( chi_sub < 0.1) cout << chi_sub;

    double Vs = U*U * chi_sub / (1 - U*chi_sub) 
        + pow(U,3)*chi_sub*chi_sub / (1 - U*U * chi_sub*chi_sub);
    return Vs;
}

double potential_scalapino_cube(Vec k1, Vec k2, double w, double T, const unordered_map<double, vector<vector<vector<double>>>> &chi_map) {
    Vec q_minus = to_IBZ_2(k1 - k2);
    Vec q_plus = to_IBZ_2(k1 + k2);
    w = round(w,6);
    //if (w!=0) cout << w << endl;

    auto chi_cube = chi_map.at(w);

    double chi_minus = calculate_chi_from_cube(chi_cube, q_minus);
    double chi_plus = calculate_chi_from_cube(chi_cube, q_plus);

    double V_minus = U*U * chi_minus / (1 - U*chi_minus) 
        + pow(U,3)*chi_minus*chi_minus / (1 - U*U * chi_minus*chi_minus);
    double V_plus = U*U * chi_plus / (1 - U*chi_plus) 
        + pow(U,3)*chi_plus*chi_plus / (1 - U*U * chi_plus*chi_plus);

    //return V_minus;
    return 0.5 * (V_minus + V_plus);
}

double potential_scalapino_triplet(Vec k1, Vec k2, double w, double T, const unordered_map<double, vector<vector<vector<double>>>> &chi_map) {
    Vec q_minus = to_IBZ_2(k1 - k2);
    Vec q_plus = to_IBZ_2(k1 + k2);

    auto chi_cube = chi_map.at(w);

    double chi_minus = calculate_chi_from_cube(chi_cube, q_minus);
    double chi_plus = calculate_chi_from_cube(chi_cube, q_plus);

    double V_minus = -pow(U,2) * chi_minus / ( 1 - pow(U*chi_minus,2));
    double V_plus = -pow(U,2) * chi_plus / ( 1 - pow(U*chi_plus,2));

    return V_minus;
    return 0.5 * ( V_minus - V_plus);
}

// Scalapino Potential Section
double get_k(double i, double n) {
    return k_max*(2.0*i/(n-1.0)-1.0);
}

double f(double E, double T) {
    if (T == 0) {
        if (E < 0) return 1;
        return 0;
    }
    return 1 / (1 + exp(E/T));
}

double ratio(Vec q, Vec k, double T, double mu, double w) {
    double e_qk = epsilon(q+k) - mu;
    double e_k = epsilon(k) - mu;
    double f_k = f(e_k, T);
    double f_qk = f(e_qk, T);

    double dE = e_qk - e_k;
    if (fabs(dE) < 0.0001 and fabs(w) < 0.0001) {
        if (T==0 or exp(e_k/T) > 1e6) {
            return e_k < 0;
        }
        double term1 = 1/T * exp(e_k/T) / pow( exp(e_k/T) + 1,2);
        //cout << "Ratio: " << term1 << endl;
        return term1;
    }
    return (f_qk - f_k) / (e_k - e_qk + w);
}

double modified_ratio(Vec q, Vec k, double T, double mu, double w, double delta) {
    double e_plus = epsilon(k + q/2) - mu;
    double e_minus = epsilon(k - q/2) - mu;
    double f_minus = f(e_minus, T);
    double f_plus = f(e_plus, T);
    double dE = e_plus - e_minus;

    if (abs(dE - w) < delta) {
        //cout << "Modified Ratio: " << 0.5 * 1/T * 1/(cosh((epsilon(k)-mu)/T) + 1) << endl;
        if (w == 0) return 0.5 * 1/T * 1/(cosh((epsilon(k)-mu)/T) + 1);
        double ep_base = e_base(k,q) - mu;
        double num = sinh(delta / T) * (cosh(ep_base / T)*cosh(w / (2*T)) + cosh(delta / T));
        double denom = (cosh((w/2 - delta)/T) + cosh(ep_base/T))*(cosh((w/2 + delta)/T) + cosh(ep_base/T));
        double avg = num / (2 * delta * denom);
        if (isnan(avg)) return 0;
        return avg;
    }
    return (f_plus - f_minus) / (w - dE);
}

double imaginary_ratio(Vec q, Vec k, double T, double mu, double w, double eta) {
    double e_plus = epsilon(k + q/2) - mu;
    double e_minus = epsilon(k - q/2) - mu;
    double f_minus = f(e_minus, T);
    double f_plus = f(e_plus, T);
    double dE = e_plus - e_minus;

    if (q.norm() < 0.0001) {
        if (w == 0) return 0.5 * 1/T * 1/(cosh((epsilon(k)-mu)/T) + 1);
        return 0;
       // double ep_base = e_base(k,q);
       // double num = sinh(delta / T) * (cosh(ep_base / T)*cosh(w / (2*T)) + cosh(delta / T));
       // double denom = (cosh((w/2 - delta)/T) + cosh(ep_base/T))*(cosh((w/2 + delta)/T) + cosh(ep_base/T));
       // double avg = num / (2 * delta * denom);
       // if (isnan(avg)) return 0;
       // return avg;
    }

    return (f_plus - f_minus) * (w - dE) / ( pow(w - dE,2) + pow(eta,2));
}


double imaginary_integration(Vec q, double T, double mu, double w, int num_points, double eta) {
    auto func = [q, T, mu, w, eta] (double kx, double ky, double kz) {
        //return modified_ratio(q, Vec(kx, ky, kz), T, mu, w, 0.0001);
        return imaginary_ratio(q, Vec(kx, ky, kz), T, mu, w, eta);
    };
    
    double x0 = -k_max + q.vals[0]/2, x1 = k_max + q.vals[0]/2;
    double y0 = -k_max + q.vals[1]/2, y1 = k_max + q.vals[1]/2;
    double z0 = -k_max + q.vals[2]/2, z1 = k_max + q.vals[2]/2;

    return trapezoidal_integration(func, x0, x1, y0, y1, z0, z1, num_points) / pow(2*k_max,dim);
}

double chi_trapezoidal(Vec q, double T, double mu, double w, int num_points) {
    double sum = 0;
    #pragma omp parallel for reduction(+:sum)
    //int num_skipped = 0;
    //ofstream file("chi_temp2.txt");
    for (int i = 0; i < num_points; i++) {
        double temp = i;
        double x = get_k(temp, num_points); 
        for (double j = 0; j < num_points; j++) {
            double y = get_k(j, num_points);
            for (double k = 0; k < num_points * (dim%2) + 1 * ((dim+1)%2); k++) {
                double z = get_k(k, num_points);
                double weight = 1.0;
                if (i == 0 or i == num_points - 1) weight /= 2.0;
                if (j == 0 or j == num_points - 1) weight /= 2.0;
                if ( (k == 0 or k == num_points - 1) and dim == 3) weight /= 2.0;

                Vec k_val(x, y, z);

                double r = ratio(q, k_val, T, mu, w);
                //cout << k_val << r << endl;
                sum += weight*r;
            }
        }
    }
    return sum / pow(num_points-1,dim); 
}

double integrate_susceptibility(Vec q, double T, double mu, double w, int num_points) {
    //return imaginary_integration(q, T, mu, w, num_points, 0.000);
    if (q.norm() < 0.0001) {
        vector<Vec> FS = tetrahedron_method(e_base_avg, q, mu);
        double sum = 0; for (auto x : FS) sum += x.area / vp(x);
        return sum / pow(2*k_max,dim);
    }
    if (q.norm() < 0.0001) q = Vec(0.01,0.01,0.01);

    double a, b; get_bounds3(q, b, a, denominator);
    vector<double> spacing; get_spacing_vec(spacing, w, a, b, num_points);
    double c = tetrahedron_sum_continuous(denominator, denominator_diff, q, spacing, w, T);
    return c / pow(2*k_max,dim);

    auto func = [q, T, mu, w] (double kx, double ky, double kz) {
        //if (epsilon(Vec(kx,ky,kz)) < mu) return 1;
        //return 0;
        return ratio(q, Vec(kx, ky, kz), T, mu, w);
        //return modified_ratio(q, Vec(kx, ky, kz), T, mu, w, 0.0001);
    };

    int base_div = 20;
    int x_divs = base_div, y_divs = base_div, z_divs = base_div;
    if (dim == 2) z_divs = 1;
    double x0 = -k_max + q.vals[0]/2, x1 = k_max + q.vals[0]/2;
    double y0 = -k_max + q.vals[1]/2, y1 = k_max + q.vals[1]/2;
    double z0 = -k_max + q.vals[2]/2, z1 = k_max + q.vals[2]/2;
    //return adaptive_trapezoidal(func, -k_max, k_max, -k_max, k_max, -k_max, k_max, x_divs, y_divs, z_divs, 0.01) / pow(2*k_max,dim);
    //return chi_trapezoidal(q, T, mu, 100);
    //return modified_integration(func, -k_max, k_max, -k_max, k_max, -k_max, k_max, 100, 0, 0, 0);
    return trapezoidal_integration(func, x0, x1, y0, y1, z0, z1, num_points) / pow(2*k_max,dim);
}

double trapezoidal_integration(auto &f, double x0, double x1, double y0, double y1, double z0, double z1, int num_points) {
    double sum = 0;
    double dx = (x1 - x0) / (num_points - 1);
    double dy = (y1 - y0) / (num_points - 1);
    double dz = (z1 - z0) / (num_points - 1);
    if (dim == 2) dz = 1;
    int counter = 0;
    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < num_points; i++) {
        double x = x0 + i*dx;
        for (int j = 0; j < num_points; j++) {
            double y = y0 + j*dy;
            for (double k = 0; k < num_points * (dim%2) + 1 * ((dim+1)%2); k++) {
                double z = z0 + k*dz;
                double w = 1.0;
                if (i == 0 or i == num_points - 1) w /= 2.0;
                if (j == 0 or j == num_points - 1) w /= 2.0;
                if ( (k == 0 or k == num_points - 1) and dim == 3) w /= 2.0;

                //if (x*x + y*y + z*z > k_max*k_max) continue;
                //cout << x << " " << y << " " << z << " " << f(x,y,z) << endl;
                sum += w * f(x,y,z);
            }
        }
    }
    return sum * dx * dy * dz;
}

double modified_integral_wrapper(Vec q, double T, double mu, double w, double delta, int num_points) {
    auto func = [q, T, mu, w, delta] (double kx, double ky, double kz) {
        return modified_ratio(q, Vec(kx, ky, kz), T, mu, w, delta);
    };
    auto s_func = [q, mu, w] (double kx, double ky, double kz) {
        Vec k(kx, ky, kz);
        return w - e_split(k,q);
    };

    double x0 = -k_max + q.vals[0]/2, x1 = k_max + q.vals[0]/2;
    double y0 = -k_max + q.vals[1]/2, y1 = k_max + q.vals[1]/2;
    double z0 = -k_max + q.vals[2]/2, z1 = k_max + q.vals[2]/2;

    return trapezoidal_integration(func, x0, x1, y0, y1, z0, z1, num_points) / pow(2*k_max,dim);
    //return modified_integration(func, -k_max, k_max, -k_max, k_max, -k_max, k_max, 100, s_func, delta, avg) / pow(2*k_max,dim);
}

double modified_integral_wrapper_1D(double a, double b, double delta, int num_points) {
    auto func = [delta] (double x) {
        return x / (x*x + delta*delta);
    };
    auto s_func = [] (double x) {
        return x;
    };
    auto avg = [] (double x) {
        return 0;
    };
    return modified_integration_1D(func, a, b, num_points, s_func, delta, avg);
}

double modified_integration(auto &f, double x0, double x1, double y0, double y1, double z0, double z1, int num_points, auto &surface, double delta, auto &avg) {
    double sum = 0;
    double dx = (x1 - x0) / (num_points - 1);
    double dy = (y1 - y0) / (num_points - 1);
    double dz = (z1 - z0) / (num_points - 1);
    if (dim == 2) dz = 1;
    #pragma omp parallel for reduction(+:sum) 
    for (int i = 0; i < num_points; i++) {
        double x = x0 + i*dx;
        for (int j = 0; j < num_points; j++) {
            double y = y0 + j*dy;
            for (double k = 0; k < num_points * (dim%2) + 1 * ((dim+1)%2); k++) {
                double z = z0 + k*dz;
                double w = 1.0;
                if (i == 0 or i == num_points - 1) w /= 2.0;
                if (j == 0 or j == num_points - 1) w /= 2.0;
                if ( (k == 0 or k == num_points - 1) and dim == 3) w /= 2.0;

                double val = f(x,y,z); 
                if (abs(surface(x,y,z)) < delta) {
                    //cout << "Surface Hit\n";
                    val = avg(x,y,z);
                    //val = 0;
                }

                sum += w * val; 
            }
        }
    }
    return sum * dx * dy * dz; 
}

double modified_integration_1D(auto &f, double x0, double x1, int num_points, auto &surface, double delta, auto &avg) {
    double sum = 0;
    double dx = (x1 - x0) / (num_points - 1);
    #pragma omp parallel for reduction(+:sum) 
    for (int i = 0; i < num_points; i++) {
        double x = x0 + i*dx;
        double w = 1.0;
        if (i == 0 or i == num_points - 1) w /= 2.0;

        sum += w * f(x) * dx; 
    }
    return sum;
}

double trap_cube(auto &f, double x0, double x1, double y0, double y1, double z0, double z1) {
    return (x1-x0)*(y1-y0)*(z1-z0) / 8 * (f(x0,y0,z0) + f(x0,y0,z1) + f(x0,y1,z0) + f(x0,y1,z1) + f(x1,y0,z0) + f(x1,y0,z1) + f(x1,y1,z0) + f(x1,y1,z1));
}

double trap_8_cubes(auto &f, double x0, double x1, double y0, double y1, double z0, double z1) {
    return trap_cube(f, x0, (x0+x1)/2, y0, (y0+y1)/2, z0, (z0+z1)/2) 
        + trap_cube(f, x0, (x0+x1)/2, y0, (y0+y1)/2, (z0+z1)/2, z1) 
        + trap_cube(f, x0, (x0+x1)/2, (y0+y1)/2, y1, z0, (z0+z1)/2) 
        + trap_cube(f, x0, (x0+x1)/2, (y0+y1)/2, y1, (z0+z1)/2, z1) 
        + trap_cube(f, (x0+x1)/2, x1, y0, (y0+y1)/2, z0, (z0+z1)/2) 
        + trap_cube(f, (x0+x1)/2, x1, y0, (y0+y1)/2, (z0+z1)/2, z1) 
        + trap_cube(f, (x0+x1)/2, x1, (y0+y1)/2, y1, z0, (z0+z1)/2) 
        + trap_cube(f, (x0+x1)/2, x1, (y0+y1)/2, y1, (z0+z1)/2, z1);
}

double adaptive_trapezoidal(auto &f, double x0, double x1, double y0, double y1, double z0, double z1, int xdivs, int ydivs, int zdivs, double error_relative) {
    double sum = 0;

    double dx = (x1 - x0) / (xdivs);
    double dy = (y1 - y0) / (ydivs);
    double dz = (z1 - z0) / (zdivs);
    if (dim == 2) dz = 1;

    for (int i = 0; i < xdivs; i++) {
        double x = x0 + i*dx;
        for (int j = 0; j < ydivs; j++) {
            double y = y0 + j*dy;
            for (double k = 0; k < zdivs; k++) {
                double z = z0 + k*dz;

                double t1 = trap_cube(f, x, x+dx, y, y+dy, z, z+dz);
                double t2 = trap_8_cubes(f, x, x+dx, y, y+dy, z, z+dz);

                if (fabs(t1 - t2) < error_relative * fabs(t2) or fabs(t1 - t2) < 0.0001) {
                    sum += t2;
                }
                else {
//                    cout << t1 << " " << t2 << " " << fabs(t1 - t2) << " " << error_relative * fabs(t2) << endl;
                    double new_zdiv = 2 * (dim % 2) + 1 * ((dim+1)%2);
                    sum += adaptive_trapezoidal(f, x, x+dx, y, y+dy, z, z+dz, 2, 2, new_zdiv, error_relative);
                }

            }
        }
    }
    return sum;
}

double iteratively_splitting_cubes(auto &f, double x0, double x1, double y0, double y1, double z0, double z1, double error_total, double error_relative) {

    double total_sum = 0;
    double dx = (x1 - x0) / 2;
    bool no_errors = false;

    for (int i = 1; not no_errors; i++) {
        total_sum = 0;
        no_errors = true;
        int iters = pow(2,i);
        //#pragma omp parallel for reduction(+:total_sum)
        for (int j = 0; j < iters; j++) {
            for (int k = 0; k < iters; k++) {
                for (int l = 0; l < iters; l++) {
                    double t1 = trap_cube(f, x0+j*dx, x0+(j+1)*dx, y0+k*dx, y0+(k+1)*dx, z0+l*dx, z0+(l+1)*dx);
                    double t2 = trap_8_cubes(f, x0+j*dx, x0+(j+1)*dx, y0+k*dx, y0+(k+1)*dx, z0+l*dx, z0+(l+1)*dx); 
                    double err = fabs(t2 - t1);
                    if ( err > fabs(error_relative*t2) and err > error_total / pow(2,3*i) )
                        no_errors = false;

                    total_sum += t2;
                }
            }
        }
        //cout << "Splits: " << i << endl;
        dx /= 2;
    }
    return total_sum;
}

vector<vector<vector<double>>> chi_cube(double T, double mu, double DOS, double w) {
    int m_z = m*(dim%2) + 3*((dim+1)%2);
    vector<vector<vector<double>>> cube(m, vector<vector<double>> (m, vector<double> (m_z)));
    unordered_map<string, double> map;
    cout << "Calculating Chi Cube...\n";
    double empty_val = -98214214;
    map[vec_to_string(Vec(0,0,0))] = DOS;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            for (int k = 0; k < m_z; k++) {
                Vec q((2*k_max*i)/(m-1), (2*k_max*j)/(m-1), (2*k_max*k)/(m_z-1));
                Vec q2 = to_IBZ_2(q);
                if (map.find(vec_to_string(q2)) == map.end())
                    map[vec_to_string(q2)] = empty_val;
            }
        }
    }
    
    cout << "Taking " << map.size() << " integrals in " << dim << " dimensions.\n";
    //#pragma omp parallel for
    for(unsigned int i = 0; i < map.size(); i++) {
        auto datIt = map.begin();
        advance(datIt, i);
        string key = datIt->first;
        map[key] = integrate_susceptibility(string_to_vec(key), T, mu, w, s_pts);
        progress_bar(1.0 * i / map.size());
    }

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            for (int k = 0; k < m_z; k++) {
                Vec q((2*k_max*i)/(m-1), (2*k_max*j)/(m-1), (2*k_max*k)/(m_z-1));
                Vec q2 = to_IBZ_2(q);
                cube[i][j][k] = map[vec_to_string(q2)];
            }
        }
    }

    cout << "\nChi Cube Created.\n";
    return cube;
}


double calculate_chi_from_cube(const vector<vector<vector<double>>> &chi_cube, Vec q) {
    Vec v = to_IBZ_2(q);
    double d = 2*k_max/(m-1);

    double x = v.vals[0], y = v.vals[1], z = v.vals[2];
    if (dim == 2) z = 0;

    int i = floor(x / d);
    int j = floor(y / d);
    int k = floor(z / d);

    double x1 = i * d; 
    double y1 = j * d; 
    double z1 = k * d; 

    double x2 = x1 + d; 
    double y2 = y1 + d; 
    double z2 = z1 + d; 

    double dx = 0, dy = 0, dz = 0, wx = 0, wy = 0, wz = 0, w0 = 0;

    // Make sure there's no issue with indexing
    //cout << q << q.vals[2] << endl;
    //int s = chi_cube.size()-1; 
    //assert( i < s and j < s and k < chi_cube[0][0].size()-1);

    double f1 = chi_cube[i][j][k], f2 = chi_cube[i+1][j][k];
    double f3 = chi_cube[i+1][j+1][k], f4 = chi_cube[i][j+1][k];
    double f5 = chi_cube[i][j][k+1], f6 = chi_cube[i+1][j][k+1];
    double f7 = chi_cube[i+1][j+1][k+1], f8 = chi_cube[i][j+1][k+1];

    if (x - x1 <= z2 - z and x - x1 >= y - y1) {// blue @ 1
        w0 = f1;
        wx = (f2 - f1) / d; dx = x - x1;
        wy = (f3 - f2) / d; dy = y - y1;
        wz = (f5 - f1) / d; dz = z - z1;
    }

    else if (y + z <= z1 + y2 and x - x1 <= y - y1) {// orange @ 1
        w0 = f1;
        wx = (f3 - f4) / d; dx = x - x1;
        wy = (f4 - f1) / d; dy = y - y1;
        wz = (f5 - f1) / d; dz = z - z1;
    }

    else if (x + z >= z1 + x2 and y + z <= z1 + y2) {// red @ 2
        w0 = f2;
        wx = (f6 - f5) / d; dx = x - x2;
        wy = (f3 - f2) / d; dy = y - y1;
        wz = (f6 - f2) / d; dz = z - z1;
    }

    else if (y + z >= z1 + y2 and x + z <= z1 + x2) {// purple @ 4
        w0 = f4;
        wx = (f3 - f4) / d; dx = x - x1;
        wy = (f8 - f5) / d; dy = y - y2;
        wz = (f8 - f4) / d; dz = z - z1;
    }

    else if (x - x1 >= y - y1 and y + z >= z1 + y2) {// teal @ 7
        w0 = f7;
        wx = (f6 - f5) / d; dx = x - x2;
        wy = (f7 - f6) / d; dy = y - y2;
        wz = (f7 - f3) / d; dz = z - z2;
    }

    else if (x - x1 <= y - y1 and y + z >= z1 + y2) {// green @ 7
        w0 = f7;
        wx = (f7 - f8) / d; dx = x - x2;
        wy = (f8 - f5) / d; dy = y - y2;
        wz = (f7 - f3) / d; dz = z - z2;
    }
    else return f1 + (f2-f1)/d*x + (f4-f1)/d*y + (f5-f1)/d*z;

    return w0 + wx*dx + wy*dy + wz*dz;
}

Vec to_IBZ_2(const Vec k) {
    Vec q = k;
    if (q.cartesian == false) q.to_cartesian();
    double x = q.vals[0], y = q.vals[1], z = q.vals[2];
    x = abs(x); y = abs(y); z = abs(z);
    if (x > M_PI) x = - (x - 2*M_PI);
    if (y > M_PI) y = - (y - 2*M_PI);
    if (z > M_PI) z = - (z - 2*M_PI);
    if (dim == 3) {
        double arr[] = {x, y, z};
        sort(arr, arr+3, greater<double>());
        auto& [a, b, c] = arr;
        Vec result(a, b, c);
        return result;
    }
    else if (dim == 2) {
        double arr[] = {x, y};
        sort(arr, arr+2, greater<double>());
        auto& [a, b] = arr;
        Vec result(a, b, z);
        return result;
    }
    else {
        cout << "Wrong Dimension\n";
        return q;
    }
}

Vec to_IBZ_spherical(const Vec k) {
    double theta = k.vals[1];
    double phi = 0;
    if (dim == 3) phi = k.vals[2];

    if (theta < 0) theta += 2*M_PI; //{cout << theta << ", "; theta += 2*M_PI; cout << theta << ", ";}
    if (theta <= M_PI/2 and theta > M_PI/4) theta = M_PI/2 - theta;
    if (theta <= 3*M_PI/4 and theta > M_PI/2) theta = theta - M_PI/2;
    if (theta <= M_PI and theta > 3*M_PI/4) theta = M_PI - theta;
    if (theta <= 5*M_PI/4 and theta > M_PI) theta = theta - M_PI;
    if (theta <= 6*M_PI/4 and theta > 5*M_PI/4) theta = 3*M_PI/2 - theta;
    if (theta <= 7*M_PI/4 and theta > 6*M_PI/4) theta = theta - 3*M_PI/2;
    if (theta <= 2*M_PI and theta > 7*M_PI/4) theta = 2*M_PI - theta;
    //if (theta < 0) cout << theta << endl;


    //if (phi <= M_PI/2 and phi > M_PI/4) {cout << "FIRST"; phi = M_PI/2 - phi;}
    //if (phi <= -M_PI/4 and phi > -M_PI/2) {cout << "THIRD" ; phi = M_PI/2 + phi;}

    Vec q(k.vals[0], theta, phi, false);
    return q;
}

