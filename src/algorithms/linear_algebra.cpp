#include <ctime>
#include <iostream>
#include <unistd.h>

#include <boost/functional/hash.hpp>
#include <math.h>
#include <openblas/lapacke.h>
#include <vector>

#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/tools/roots.hpp>

#include "../config/load/cpp_config.hpp"
#include "../objects/vec.hpp"
#include "../objects/eigenvec.hpp"
#include "../objects/matrix.hpp"

using namespace std;

/* Power iteration algorithm
 * Not optimized, but so fast it doesn't matter
 * Returns eigenvalue and eigenvector
 * NOTE: Technically returns multiple (up to 2) eigenvalues and eigenvectors
 * If the top 2 eigenvalues are positive and nonzero it will return 2 solutions instead of 1
 */
vector<Eigenvector> power_iteration(Matrix &A, float error) {
    printf("power iteration called\n");
    vector<Eigenvector> vals;
    for (int eig_num = 0; eig_num < 5; eig_num++) {
        Eigenvector x(A.size, true);
        printf("x created\n");
        float diff_mag = 1;
        float rayleigh = dot(x, A * x) / dot(x, x);
        float sum = 0;

        int iterations = 0;
        printv("Eig Num: %d\n", eig_num);
        for (int i = 0; diff_mag > error; i++) {
            Eigenvector x_new = A * x;
            x_new.normalize();

            // Account for phase shift
            float cos_diff = dot(x, x_new) / (x.norm() * x_new.norm());
            diff_mag = abs(cos_diff) - 1.0;
            //Eigenvector diff_vec = x - x_new;
            //Eigenvector diff_neg_vec = x + x_new;
            //diff_mag = diff_vec.norm();
            //if (diff_mag > diff_neg_vec.norm())
            //    diff_mag = diff_neg_vec.norm();

            // Update iterated eigenvector
            // if ( i%100 == 0) cout << x.transpose()*A*x << endl;
            x = x_new;
            iterations = i;
            // cout << i << " " << 1000 << " " << i%1000 << endl;
            if (i % 1000 == 0)
                printv("Iteration: %d %f\n", i, dot(x, A * x));
            if (i % 100 == 0 and i > 50000 and dot(x, A * x) < 0)
                break;
        }
        printv("Iterations: %d\n", iterations);
        rayleigh = dot(x, A * x) / dot(x, x);
        x.eigenvalue = rayleigh + sum;

        vals.push_back(x);
        if (x.eigenvalue > 0)
            return vals;

        sum += x.eigenvalue;
        Matrix identity(A.size);
        A = A - identity * x.eigenvalue;
        Eigenvector temp(A.size);
        x = temp;
        diff_mag = 1;
    }
    return vals;
}

// This function is designed to get the eigenvalue accurately, not the eigenvector
Eigenvector power_iteration(Matrix &A, Eigenvector initial) {
    Eigenvector x = initial;
    //Eigenvector x(A.size, true);
    float diff_percent = 1;
    vector<float> diffs(10, 1.0f);
    bool all_low = false;
    for (int i = 0; !all_low; i++) {
        Eigenvector x_new = A * x;
        x_new.normalize();
        x_new.eigenvalue = dot(x_new, A * x_new) / dot(x_new, x_new);
        float cos_diff = dot(x, x_new) / (x.norm() * x_new.norm());
        //diff_percent = abs(abs(cos_diff) - 1.0);
        diff_percent = abs(x.eigenvalue - x_new.eigenvalue) / (abs(x.eigenvalue) + 1e-6);
        diffs[i % 10] = diff_percent;
        if (i % 10 == 0) {
            printv("Iteration: %d, Eigenvalue: %f, Diff Percent: %f\n", i, x_new.eigenvalue, diff_percent);
            all_low = true;
            for (int j = 0; j < 10; j++)
                if (diffs[j] > 1e-4)
                    all_low = false;
        }
        x = x_new;
    }
    if (x.eigenvalue < 0) {
        A -= Matrix(A.size) * x.eigenvalue;
        return power_iteration(A, initial);
    }
    return x;
}

void lapack_diagonalization(Matrix &A, Eigenvector *eigenvectors) {
    int N = A.size;

    // Allocate memory for eigenvalues and eigenvectors
    float *val_r = new float[N];    // Real parts of eigenvalues
    float *val_i = new float[N];    // Imaginary parts of eigenvalues (for
                                    // non-symmetric matrices)
    float *vecs = new float[N * N]; // Eigenvectors

    int info;

    // LAPACKE_sgeev parameters:
    // 'N' -> No computation of left eigenvectors
    // 'V' -> Compute right eigenvectors
    info = LAPACKE_sgeev(LAPACK_ROW_MAJOR, 'N', 'V', N, A.vals, N, val_r, val_i,
                         NULL, N, vecs, N);

    if (info > 0) {
        std::cerr << "The algorithm failed to compute eigenvalues."
                  << std::endl;
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
    const int il =
        N + 1 -
        num_eigenvalues_to_save; // Index of the first eigenvalue to be found
    const int iu = N;            // Index of the last eigenvalue to be found
    const char jobz = 'V';       // Compute eigenvalues and eigenvectors
    const char range = 'I';      // Compute a subset of eigenvalues
    const char uplo = 'L';       // Lower triangle of A is stored
    const float abstol = 0.0f;   // Use default tolerance

    float *w = new float[N];     // Array for eigenvalues
    float *z = new float[N * N]; // Array for eigenvectors
    int isuppz[2 * N];           // Support array for eigenvectors

    int m; // Total number of eigenvalues found
    int info;

    int count = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (A.vals[i * N + j] != 0) {
                count++;
            }
        }
    }
    printf("Diagonalizing Matrix\n");
    info = LAPACKE_ssyevr(LAPACK_ROW_MAJOR, jobz, range, uplo, N, A.vals, lda,
                          0.0f, 0.0f, il, iu, abstol, &m, w, z, lda, isuppz);

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
            if (z[i * N + j] != 0) {
                count++;
            }
        }
    }

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

    delete[] w;
    delete[] z;
}
