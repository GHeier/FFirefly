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
#include <memory>
#include <lapacke.h>
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

Eigenvector power_iteration(Matrix &A) {
    Eigenvector x(A.size, true);
    float diff_percent = 1;
    for (int i = 0; diff_percent < 0.01; i++) {
        Eigenvector x_new = A * x;
        x_new.normalize();
        x_new.eigenvalue = dot(x_new, A * x_new) / dot(x_new, x_new);
        diff_percent = fabs((x_new.eigenvalue - x.eigenvalue) / x.eigenvalue);
        printf("Iteration: %d, Eigenvalue: %f, Diff Percent: %f\n", i, x_new.eigenvalue, diff_percent);
        x = x_new;
    }
    if (x.eigenvalue < 0) {
        A -= Matrix(A.size)*x.eigenvalue;
        return power_iteration(A);
    }
    return x;
}

void lapack_diagonalization(Matrix &A, Eigenvector *eigenvectors) {
    int N = A.size;
    
    // Allocate memory for eigenvalues and eigenvectors
    float *val_r = new float[N]; // Real parts of eigenvalues
    float *val_i = new float[N]; // Imaginary parts of eigenvalues (for non-symmetric matrices)
    float *vecs = new float[N * N]; // Eigenvectors

    int info;
    
    // LAPACKE_sgeev parameters:
    // 'N' -> No computation of left eigenvectors
    // 'V' -> Compute right eigenvectors
    info = LAPACKE_sgeev(LAPACK_ROW_MAJOR, 'N', 'V', N, A.vals, N, val_r, val_i, NULL, N, vecs, N);
    
    if (info > 0) {
        std::cerr << "The algorithm failed to compute eigenvalues." << std::endl;
        delete[] val_r;
        delete[] val_i;
        delete[] vecs;
        return;
    }
    
    // Convert to Eigenvector format
    for (int i = 0; i < num_eigenvalues_to_save; i++) {
        eigenvectors[i] = Eigenvector(N);
        eigenvectors[i].eigenvalue = val_r[i];
        for (int j = 0; j < N; j++) {
            eigenvectors[i][j] = vecs[i * N + j];
        }
    }

    // Clean up
    delete[] val_r;
    delete[] val_i;
    delete[] vecs;
}

void lapack_hermitian_diagonalization(Matrix &A, Eigenvector *eigenvectors) {
    const int N = A.size; // Dimension of the matrix
    const int lda = N;
    const int il = N + 1 - num_eigenvalues_to_save; // Index of the first eigenvalue to be found
    const int iu = N; // Index of the last eigenvalue to be found
    const char jobz = 'V'; // Compute eigenvalues and eigenvectors
    const char range = 'I'; // Compute a subset of eigenvalues
    const char uplo = 'L'; // Lower triangle of A is stored
    const float abstol = 0.0f; // Use default tolerance

    float *w = new float[N]; // Array for eigenvalues
    float *z = new float[N * N]; // Array for eigenvectors
    int isuppz[2 * N]; // Support array for eigenvectors

    int m; // Total number of eigenvalues found
    int info;

    int count = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (A.vals[i*N+j] != 0) {
                count++;
            }
        }
    }
    printf("Nonzero elements: %d\n", count);

    printf("Diagonalizing Matrix\n");
    info = LAPACKE_ssyevr(LAPACK_ROW_MAJOR, jobz, range, uplo, N,
                          A.vals, lda, 0.0f, 0.0f, il, iu, abstol, &m, w,
                          z, lda, isuppz);
    
    if (info != 0) {
        std::cerr << "Error: LAPACKE_ssyevr returned " << info << "\n";
        delete[] w;
        delete[] z;
        return;
    }
    
    printf("Diagonalization Successful\n");
    count = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (z[i*N+j] != 0) {
                count++;
            }
        }
    }
    printf("Nonzero elements: %d\n", count);

    count = 0;
    // Convert to Eigenvector format
    for (int i = 0; i < num_eigenvalues_to_save; ++i) {
        eigenvectors[i] = Eigenvector(N);
        eigenvectors[i].eigenvalue = w[i];
        for (int j = 0; j < N; ++j) {
            eigenvectors[i][j] = z[j * N + i];
            if (z[j * N + i] != 0) {
                count++;
            }
        }
    }
    printf("Nonzero elements: %d\n", count);

    delete[] w; 
    delete[] z;
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
        //#pragma omp parallel for
        for (int j = 0; j < P.size; j++) {
            Vec k2 = k[j];
            P(i,j) = (float)(-pow(k1.area/vp(k1),0.5) * V(k1, k2, 0, T, chi_cube2) * pow(k2.area/vp(k2),0.5));
            if (isnan(P(i,j)) == true) {  // debugging statement
		    cout << "P matrix issue. P(i, j) = " << P(i,j) << std::endl;
		    exit(0);
	    }
	    //assert(isnan(P(i,j)) == false); //ORIGINAL WORK
        }
        progress_bar(1.0 * i / P.size);
    }
    P *= (2 / pow(2*M_PI, dim));
    cout << "\nP Matrix Created\n";
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
    for (unsigned int i = 0; i < num_eigenvalues_to_save; i++) {
        for (unsigned int j = 0; j < vectors[i].size; j++) {
            Vec k = FS[j];
            vectors[i][j] /= pow(k.area/vp(k),0.5);
        }
    }
}

void freq_vector_to_wave(vector<vector<Vec>> &freq_FS, Eigenvector *vectors) {
    int size = matrix_size_from_freq_FS(freq_FS);
    for (unsigned int x = 0; x < num_eigenvalues_to_save; x++) {
        int ind = 0;
        for (unsigned int i = 0; i < freq_FS.size(); i++) {
            for (unsigned int j = 0; j < freq_FS[i].size(); j++) {
                Vec k = freq_FS[i][j];
                vectors[x][ind] /= pow(k.area/vp(k),0.5);
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
