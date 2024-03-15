#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

#include <vector>
#include <algorithm>

#include <Eigen/Dense>
//#include <lambda_lanczos/lambda_lanczos.hpp>

#include "cfg.h"
#include "fermi_surface.h"
#include "calculations.h"
#include "potential.h"
#include "save_data.h"
#include "vec.h"
#include "utilities.h"

using namespace Eigen;
using std::string;

int main() {
/* 
 * ========================================================================================
 * ======================== FERMI SURFACE CREATION AND FILE SAVING ========================
   ========================================================================================
 */
    cout << "Calculating Fermi Surface..." << endl;
    vector<Vec> FS = tetrahedron_method(mu);
    cout << "Number of points along Fermi Surface: " << FS.size() << endl;
    save_FS(FS);
    double DOS = get_DOS(FS);
    cout << "Density of States: " << DOS << endl;
/* 
 * ========================================================================================
 * =========================== CRITICAL TEMPERATURE CALCULATION  ==========================
   ========================================================================================
 */
    double T = 0.25;
    cout << setprecision(10);
    //cout << coupling_calc(FS, T) << endl;
    //T = 0.065;
    //T = get_Tc(FS);
    printf("Temperature: %f \n", T);

/* 
 * ========================================================================================
 * ========================== MATRIX CREATION AND DIAGONALIZATION  ========================
   ========================================================================================
 */
    vector<vector<vector<double>>> cube;
    if (potential_name != "test") cube = chi_cube(T, mu, DOS);

    MatrixXd P = create_P(FS, T, cube);
    double f = f_singlet_integral(T);
    cout << "F-integral value: " << f << endl;

    cout << "Finding Eigenspace..." << endl;
    // Solving every vector using Eigen method
    EigenSolver<MatrixXd> s(P);
    
    VectorXcd vals = s.eigenvalues();// * f / FS.size();
    //VectorXcd vals = s.eigenvalues() * f / FS.size();
    EigenSolver<MatrixXd>::EigenvectorsType vecs;
    vecs = s.eigenvectors(); 
    
    cout << "Eigenspace Found\n";

    // Testing to confirm Eigen didn't mess up the first vector at least
    VectorXd first_vec = vecs.col(0).real();
    if ( ( P*first_vec 
                - vals(0).real() * first_vec ).norm() > 0.00001 ) {
        cout << "Matrix Decomposition Failed\n";
    }
    double mag = first_vec.transpose() * first_vec;
    if ( fabs(mag - 1.0) > 0.01 ) cout << "Eigenvector not normalized\n";
    //vals = vals * f;

    cout << "Saving Potential and Susceptibility Functions\n";
    save_potential_vs_q(FS, P, "potential.dat");
    if (cube.size() != 0 ) save_chi_vs_q(cube, FS, "chi.dat");

/* 
 * ========================================================================================
 * =========================== EIGENVALUE/VECTOR SORTING AND SAVING  ======================
   ========================================================================================
 */
    // Sort solutions with highest eigenvalue/eigenvector pair first
    std::vector<EigAndVec> solutions;
    solutions = combine_eigs_and_vecs(vals.real(), vecs.real());
    sort(solutions.rbegin(), solutions.rend());
    vector_to_wave(FS, solutions);
    //for (int i = 0; i < 6; i++) cout << solutions[0].vec(i) << " "; cout << endl;
    //vector_to_wave(FS, solutions);
    //for (int i = 0; i < 6; i++) cout << solutions[0].vec(i) << " "; cout << endl;
    
    // Defining file name based on cfg (config)
    std::ostringstream out;
    out.precision(1);
    out << std::fixed << "../data/" + potential_name << dim << "D" 
        << "_mu=" << mu << "_U=" << U << "_wD=" << w_D 
        << "_n=" << n << ".dat";
    string file_name = std::move(out).str();

    // Save file in cartesian coordinates for the sake of plotting easier
    save(file_name, T, FS, solutions);
    /*
     Using Spectre/LambdaLanczos algorithm to find only first few eigenvectors, (not that much faster)
        auto mv_mul = [&](const vector<double>& in, vector<double>& out) {
            auto eigen_in = Eigen::Map<const Eigen::VectorXd>(&in[0], in.size());
            auto eigen_out = Eigen::Map<Eigen::VectorXd>(&out[0], out.size());
            eigen_out.noalias() += V * eigen_in; // Efficient version
        };

    LambdaLanczos<double> engine(mv_mul, size, true, 3);
    vector<double> eigenvalues;
    vector<vector<double>> eigenvectors;
    engine.run(eigenvalues, eigenvectors);
    */


    return 0;
}
