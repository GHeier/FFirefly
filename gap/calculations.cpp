/*
 * @file calculations.cpp
 *
 * @brief Contains the main functions for calculating and diagonalizing the interaction matrix, 
 * critical temperature, and coupling constant
 *
 * @author Griffin Heier
 */

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

//#include <lambda_lanczos/lambda_lanczos.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/quadrature/gauss.hpp>

#include "calculations.h"
#include "cfg.h"
#include "fermi_surface.h"
#include "potential.h"
#include "vec.h"
#include "matrix.hpp"
#include "eigenvec.hpp"
#include "utilities.h"
#include "frequency_inclusion.hpp"

using std::cout;
using std::endl;
using std::sort;
using std::vector;
using std::unordered_map;
//using lambda_lanczos::LambdaLanczos;


namespace std {
    template<> struct hash<::Vec> : boost::hash<::Vec> {};
}

using vec_map = std::unordered_map<Vec, float>;

/* Power iteration algorithm
 * Not optimized, but so fast it doesn't matter
 * Returns eigenvalue and eigenvector
 * NOTE: Technically returns multiple (up to 2) eigenvalues and eigenvectors
 * This is on the off chance that all eigenvalues are positive and nonzero it will 
 * return 2 solutions instead of 1
*/ 
vector<Eigenvector> power_iteration(Matrix A, float error) {
    vector<Eigenvector> vals;
    for (int eig_num = 0; eig_num < 5; eig_num++) {
        Eigenvector x(A.size, true);
        float diff_mag = 1;
        float rayleigh = dot(x, A * x) / dot(x, x);
        float sum = 0;

        int iterations = 0;
        cout << "Eig Num: " << eig_num << endl;
        string filename = "matrix" + std::to_string(eig_num) + ".txt";
        ofstream file(filename);
        for (int i = 0; i < A.size; i++) {
            for (int j = 0; j < A.size; j++) {
                file << A(i,j) << " ";
            }
            file << endl;
        }
        for (int i = 0; diff_mag > error; i++) {
            Eigenvector x_new = A * x;
            x_new.normalize();

            // Account for phase shift
            Eigenvector diff_vec = x - x_new;
            Eigenvector diff_neg_vec = x + x_new;
            diff_mag = diff_vec.norm();
            if (diff_mag > diff_neg_vec.norm()) diff_mag = diff_neg_vec.norm();

            // Update iterated eigenvector
            //if ( i%100 == 0) cout << x.transpose()*A*x << endl;
            x = x_new;
            iterations = i;
            //cout << i << " " << 1000 << " " << i%1000 << endl;
            if (i%1000 == 0) cout << "Iteration: " << i << " " << dot(x, A * x) << endl;
            if (i%100 == 0 and i > 50000 and dot(x, A * x) < 0) break;
        }
        cout << "Iterations: " << iterations << endl;
        rayleigh = dot(x, A * x) / dot(x, x);
        x.eigenvalue = rayleigh + sum;

        vals.push_back(x);
        if (x.eigenvalue > 0) return vals;

        sum += x.eigenvalue;
        Matrix identity(A.size);
        A = A - identity*x.eigenvalue;
        Eigenvector temp(A.size);
        x = temp;
        diff_mag = 1; 
    }
    return vals;
}

Eigenvector* lapack_diagonalization(Matrix &A) {
    // All variable definitions for LAPACK
    printf("Diagonalization Step 1\n");
    int N = A.size;
    float val_r[N], val_i[N];
    float *vecs = new float[N*N];
    float work_test[1];
    int LWORK = -1;
    int info;
    char jobvl = 'N';
    char jobvr = 'V';
    printf("Diagonalization Step 2\n");
    sgeev_(&jobvl, &jobvr, &N, A.vals, &N, val_r, val_i, NULL, &N, vecs, &N, work_test, &LWORK, &info);
    printf("Diagonalization Step 3\n");
    LWORK = static_cast<int>(work_test[0]);
    float *work = new float[LWORK];
    printf("Diagonalization Step 4\n");
    sgeev_(&jobvl, &jobvr, &N, A.vals, &N, val_r, val_i, NULL, &N, vecs, &N, work, &LWORK, &info);
    printf("Diagonalization Step 5\n");
    // Convert to Eigenvec format
    Eigenvector *eigenvectors = new Eigenvector[N];
    for (int i = 0; i < N; i++) {
        Eigenvector temp(N);
        for (int j = 0; j < N; j++) {
            temp[j] = vecs[i*N+j];
        }
        temp.eigenvalue = val_r[i];
    }
    delete [] vecs;
    delete [] work;
    return eigenvectors;
}

/*
 * f_singlet is the part of the linearized BCS gap equation:
 *  epsilon * Delta_k = sum { V_kk' Delta_k' tanh(E_k / 2*T) / (2*E_k) }
 *  epsilon * Delta_k = sum { V_kk' Delta_k' f_singlet(E_k, T) }
 *  f_singlet(E_k, T) = tanh(E_k / 2*T) / (2*E_k)
 */
float f_singlet(float x, float T) {
    if (abs(x) < 0.00001) return 1/(4*T);
    return tanh(x/(2*T))/(2*x);
}

// Integral of f_singlet over E_k from -wD to wD
// wD is the debye frequency
float f_singlet_integral(float T) {
    auto f = [T](float x) {return f_singlet(x,T);};
    float integral = boost::math::quadrature::gauss<float, 7>::integrate(f, 0, wc);
    return 2*integral;
}

// Create V matrix
// Picks the potential based on the global variable "potential_name"
void create_P(Matrix &P, vector<Vec> &k, float T, const unordered_map<float, vector<vector<vector<float>>>> &chi_cube2) {
    cout << "Creating P Matrix\n";
    int a = 0;
    for (int i = 0; i < P.size; i++) {
        Vec k1 = k[i];
        #pragma omp parallel for
        for (int j = 0; j < P.size; j++) {
            Vec k2 = k[j];
            P(i,j) = (float)(-pow(k1.area/vp(k1),0.5) * V(k1, k2, 0, T, chi_cube2) * pow(k2.area/vp(k2),0.5));
            //assert(isnan(P(i,j)) == false);
        }
        progress_bar(1.0 * i / P.size);
    }
    cout << "\nP Matrix Created\n";
    P = P * (2 / pow(2*M_PI, dim));
}

// Returns the highest eigenvalue-1 of a given matrix V at temperature T
// This is the function used for root finding by get_Tc
// By finding the root of eig-1, we find the temperature where eig=1
float f(vector<Vec> k, float T) {
    cout << "\nTemperature point: " << T << endl;
    float DOS = 0; for (auto k1 : k) DOS += k1.area;
    DOS /= pow(2*M_PI, dim);
    auto cube_map = chi_cube_freq(T, mu);
    //auto cube = chi_cube(T, mu, DOS, 0);
    Matrix P(k.size());
    create_P(P, k, T, cube_map);
    float f_integrated = f_singlet_integral(T);

    vector<Eigenvector> answers = power_iteration(P, 0.0001);
    float eig = answers[answers.size() - 1].eigenvalue;
    eig *= f_integrated;

    cout << "Calculated Eigenvalue: " << eig << endl;
    return eig-1;
}

// Returns the temperature where eig=1
// Uses the f() function above to achieve that, just finds the root of the function
float get_Tc(vector<Vec> k) {
    float lower = 0.005;
    float upper = 1;

//    cout << "Determining if Tc exists...\n";
//    float max_eig = f(k, lower);
//    cout << "Maximum eigenvalue is: " << max_eig + 1 << endl;
//    assert(max_eig <= 0); 
//    cout << "Tc exists. Calculating exact Critical Temperature...\n";
    
    auto x = boost::math::tools::bisect(
            [k](float T){ return f(k,T); },
            lower,
            upper,
            [=](float lower, float upper){return upper-lower < 0.0001;}
    );

    //cout << "Lower: " << x.first << " Upper: " << x.second << endl;
    float root = (x.second + x.first) / 2;
    return root;
}

// Un-shifting the area-shifted eigenvectors in order to find wavefunction
void vector_to_wave(vector<Vec> &FS, Eigenvector *vectors) {
    for (unsigned int i = 0; i < FS.size(); i++) {
        Eigenvector temp(FS.size());
        for (unsigned int j = 0; j < FS.size(); j++) {
            Vec k = FS[j];
            temp[j] = vectors[i].eigenvector[j] / pow(k.area/vp(k),0.5);
        }
        vectors[i] = temp;
    }
}

void freq_vector_to_wave(vector<vector<Vec>> &freq_FS, Eigenvector *vectors) {
    int size = 0;
    for (auto x : freq_FS) {
        size += x.size();
    }
    for (unsigned int x = 0; x < size; x++) {
        int ind = 0;
        for (unsigned int i = 0; i < freq_FS.size(); i++) {
            for (unsigned int j = 0; j < freq_FS[i].size(); j++) {
                Vec k = freq_FS[i][j];
                vectors[x].eigenvector[ind] /= pow(k.area/vp(k),0.5);
                ind++;
            }
        }
    }
}

float get_DOS(vector<Vec> &FS) {
    float sum = 0;
    for (auto k : FS) 
        sum += k.area / vp(k);
    return sum / pow(2*M_PI, dim);
}

float coupling_calc(vector<Vec> &FS, float T) {
    cout << "Calculating Coupling Constant...\n";
    int size = FS.size();
    float DOS = get_DOS(FS);
    auto cube_map = chi_cube_freq(T, mu);
    //auto cube = chi_cube(T, mu, DOS, 0);
    float f_integrated = f_singlet_integral(T);

    float lambda = 0, normalization = 0;
    auto wave = [](Vec k) {
        Vec q = k; if (k.cartesian == false) q.to_cartesian();
        return cos(q.vals[1]) - cos(q.vals[0]);
    };
    for (int i = 0; i < size; i++) {
        Vec k1 = FS[i];
        for (int j = 0; j < size; j++) {
            Vec k2 = FS[j];
            //lambda += - k1.area * k2.area / vp(k1) * wave(k1) * V(k1, k2, T, cube) / vp(k2) * wave(k2);
            lambda += V(k1, k2, 0, T, cube_map)*k1.area*k2.area;
        }
        normalization += pow(wave(k1),2) * k1.area / vp(k1);
    }
    cout << "Normalization: " << normalization << endl;
    cout << "Lambda: " << lambda << endl;
    return lambda / normalization * (2 / pow(2*M_PI, 3)) ;
}
