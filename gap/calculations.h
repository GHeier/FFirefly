#pragma once
#ifndef CALCULATIONS_H_
#define CALCULATIONS_H_

#include <vector>
#include <Eigen/Dense>

#include "vec.h"

using namespace Eigen;
using std::vector;

struct EigAndVec {
    double eig;
    VectorXd vec;
};

struct EigsAndVecs {
    double eig;
    EigenSolver<MatrixXd>::EigenvectorsType vec;
};


vector<EigAndVec> power_iteration(MatrixXd A, double error);
double f_singlet(double x, double T);
double f_singlet_integral(double T);
MatrixXd create_P(vector<Vec> k, double T, const vector<vector<vector<double>>> &chi_cube2);
double f(vector<Vec> k, double T);
double get_Tc(vector<Vec> k);
bool operator<(const EigAndVec& left, const EigAndVec& right);
vector<EigAndVec> combine_eigs_and_vecs(VectorXd eigenvalues, EigenSolver<MatrixXd>::EigenvectorsType eigenvectors);
void vector_to_wave(vector<Vec> &FS, vector<EigAndVec> &vectors);

#endif
