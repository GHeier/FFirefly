#pragma once
#ifndef CALCULATIONS_H_
#define CALCULATIONS_H_

#include <vector>

#include "cfg.h"
#include "vec.h"
#include "matrix.hpp"
#include "eigenvec.hpp"

using std::vector;

extern "C" {
    void dgeev_( char* jobvl, char* jobvr, int* n, double* a,
                int* lda, double* wr, double* wi, double* vl, int* ldvl,
                double* vr, int* ldvr, double* work, int* lwork, int* info );
}

vector<Eigenvector> power_iteration(Matrix A, double error);
vector<Eigenvector> lapack_diagonalization(Matrix A);
double f_singlet(double x, double T);
double f_singlet_integral(double T);
void create_P(Matrix &P, vector<Vec> &k, double T, const unordered_map<double, vector<vector<vector<double>>>> &chi_cube2);
double f(vector<Vec> k, double T);
double get_Tc(vector<Vec> k);
void vector_to_wave(vector<Vec> &FS, vector<Eigenvector> &vectors);
void freq_vector_to_wave(vector<vector<Vec>> &freq_FS, vector<Eigenvector> &vectors);
double get_DOS(vector<Vec> &FS);
double coupling_calc(vector<Vec> &FS, double T);

#endif
