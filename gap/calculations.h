#pragma once
#ifndef CALCULATIONS_H_
#define CALCULATIONS_H_

#include <vector>

#include "cfg.h"
#include "vec.h"
#include "matrix.hpp"
#include "eigenvec.hpp"

using std::vector;

<<<<<<< HEAD
// LAPACK functions
extern "C" {
    void dgeev_( char* jobvl, char* jobvr, int* n, float* a,
                int* lda, float* wr, float* wi, float* vl, int* ldvl,
                float* vr, int* ldvr, float* work, int* lwork, int* info );
    void sgeev_(char* jobvl, char* jobvr, int* n, float* a,
                int* lda, float* wr, float* wi, float* vl, int* ldvl,
                float* vr, int* ldvr, float* work, int* lwork, int* info);
}

/**
 * @brief Calculates up to the first two leading eigenvectors of a matrix using the power 
 * iteration method.
 *
 * @param A The matrix to calculate the eigenvectors of.
 * @param error The error tolerance for the power iteration method.
 *
 * @return A vector of eigenvectors.
 */

vector<Eigenvector> power_iteration(Matrix A, float error);

/**
 * @brief Calculates the eigenvectors of a matrix using the LAPACK library.
 *
 * @param A The matrix to calculate the eigenvectors of.
 *
 * @return A vector of eigenvectors.
 */
Eigenvector* lapack_diagonalization(Matrix &A);

/**
 * @brief The tanh(beta*e)/(2e) function.
 *
 * @param x The energy.
 * @param T the temperature.
 *
 * @return The value of the function.
 */
float f_singlet(float x, float T);

/**
 * @brief The integral of de*tanh(beta*e)/(2e) function.
 *
 * @param T the temperature.
 *
 * @return The value of the function.
 */
float f_singlet_integral(float T);

/**
 * @brief Creates the P matrix, whose eigenvalues are the same as the gap equation
 *
 * @param P The matrix to create.
 * @param k The k points.
 * @param T The temperature.
 * @param chi_cube2 The chi matrix.
 */
void create_P(Matrix &P, vector<Vec> &k, float T, const unordered_map<float, vector<vector<vector<float>>>> &chi_cube2);

/**
 * @brief Function to be used when calculating get_Tc
 *
 * @param k The Fermi Surface points.
 * @param T The temperature.
 */
float f(vector<Vec> k, float T);

/**
 * @brief Calculates the critical temperature of a superconductor.
 *
 * Does so using the power iteration, since this calculation is done prior to finding the eigenvectors, and only needs the largest eigenvalue to find Tc
 *
 * @param k The Fermi Surface points.
 *
 * @return The critical temperature.
 */
float get_Tc(vector<Vec> k);

/**
 * @brief Modifies the eigenvectors to become the gap function.
 *
 * The calculations here transform the Matrix P to be Hermitian. This means its eigenvalues are 
 * not the gap function directly, and this function transforms it into the gap function.
 *
 * @param FS The Fermi Surface points.
 * @param vectors The eigenvectors to modify.
 */
void vector_to_wave(vector<Vec> &FS, Eigenvector *vectors);

/**
 * @brief Modifies the eigenvectors to become the gap function, for the case of defining surfaces
 * beyond the Fermi Surface.
 *
 * The calculations here transform the Matrix P to be Hermitian. This means its eigenvalues are 
 * not the gap function directly, and this function transforms it into the gap function.
 *
 * @param freq_FS The Fermi Surface points.
 * @param vectors The eigenvectors to modify.
 */
void freq_vector_to_wave(vector<vector<Vec>> &freq_FS, Eigenvector *vectors);

/**
 * @brief Calculates the density of states.
 *
 * @param FS The Fermi Surface points.
 *
 * @return The density of states.
 */
float get_DOS(vector<Vec> &FS);

/**
 * @brief Calculates the coupling constant for a d-wave
 *
 * @param FS The Fermi Surface points.
 * @param T The temperature.
 *
 * @return The coupling constant.
 */
float coupling_calc(vector<Vec> &FS, float T);
=======
vector<Eigenvector> power_iteration(Matrix &A, double error);
double f_singlet(double x, double T);
double f_singlet_integral(double T);
void create_P(Matrix &P, vector<Vec> &k, double T, const unordered_map<double, vector<vector<vector<double>>>> &chi_cube2);
double f(vector<Vec> k, double T);
double get_Tc(vector<Vec> k);
void vector_to_wave(vector<Vec> &FS, vector<Eigenvector> &vectors);
double get_DOS(vector<Vec> &FS);
double coupling_calc(vector<Vec> &FS, double T);
>>>>>>> origin/main

#endif
