/**
 * @file main.cpp
 *
 * @brief Main file for the program
 *
 * @details This file finds the Fermi Surface(s), calculates the critical temperature, finds the 
 * pairing symmetry, and saves the Gap functions to a file.
 *
 * @author Griffin Heier
 */

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
#include "frequency_inclusion.hpp"
#include "calculations.h"
#include "potential.h"
#include "save_data.h"
#include "vec.h"
#include "matrix.hpp"
#include "eigenvec.hpp"
#include "utilities.h"

using std::string;

int main() {
    // Sets the number of threads used in parallelization to one less than the maximum
    // This allows for the main thread to be used for other tasks
    int num_procs = omp_get_num_procs();
    omp_set_num_threads(num_procs - 1);






/* 
 * ========================================================================================
 * ======================== FERMI SURFACE CREATION AND FILE SAVING ========================
   ========================================================================================
 */
    cout << "Calculating Fermi Surface..." << endl;

    vector<vector<Vec>> freq_FS;
    freq_FS = freq_tetrahedron_method(mu);
    vector<Vec> FS = freq_FS[(l+1)/2 - 1];


    cout << "Number of points along Fermi Surface: " << FS.size() << endl;
    save_FS(FS);
    float DOS = get_DOS(FS);
    cout << "Density of States: " << DOS << endl;











/* 
 * ========================================================================================
 * =========================== CRITICAL TEMPERATURE CALCULATION  ==========================
   ========================================================================================
 */
    float T = 0.25;
    cout << setprecision(10);
    //cout << coupling_calc(FS, T) << endl;
    //T = 0.065;
    //T = get_Tc(FS);
    printf("Temperature: %.5f \n", T);












/* 
 * ========================================================================================
 * ========================== MATRIX CREATION AND DIAGONALIZATION  ========================
   ========================================================================================
 */
<<<<<<< HEAD
    unordered_map <float, vector<vector<vector<float>>>> cube_freq_map;
    // Calculates the susceptibility matrix if it's going to be used in the potential
    // Otherwise it's passed as empty
    if (potential_name.find("scalapino") != string::npos) 
        cube_freq_map = chi_cube_freq(T, mu);

    // This calculates total size of the matrix
    int size = matrix_size_from_freq_FS(freq_FS);

    //Matrix Pf2(size); 
    //create_P_freq(Pf2, freq_FS, T, cube_freq_map);
    Matrix P(FS.size());
    create_P(P, FS, T, cube_freq_map);

    float f = f_singlet_integral(T);
    cout << "F-integral value: " << f << endl;

    cout << "Finding Eigenspace..." << endl;
    Eigenvector *lapack_solutions = lapack_diagonalization(P);
=======
    unordered_map <double, vector<vector<vector<double>>>> cube_freq_map;
//    vector<vector<vector<double>>> cube;
//    if (potential_name != "test") cube = chi_cube(T, mu, DOS, 0);
    if (potential_name != "test" or potential_name != "const") cube_freq_map = chi_cube_freq(T, mu, DOS);
    return 0;

    int size = 0;
    for (int i = 0; i < freq_FS.size(); i++) {
        size += freq_FS[i].size();
    }
    Matrix Pf2(size); 
    create_P_freq(Pf2, freq_FS, T, cube_freq_map);
    Matrix P(FS.size());
    create_P(P, FS, T, cube_freq_map);
    double f = f_singlet_integral(T);
    cout << "F-integral value: " << f << endl;

    cout << "Finding Eigenspace..." << endl;
    vector<Eigenvector> answers = power_iteration(P, 0.001);
//    vector<EigAndVec> answersf1 = power_iteration(Pf1, 0.001);
    vector<Eigenvector> answersf2 = power_iteration(Pf2, 0.001);
    double eig = answers[answers.size() - 1].eigenvalue;
//    double eigf1 = answersf1[answersf1.size() - 1].eig;
    double eigf2 = answersf2[answersf2.size() - 1].eigenvalue;
    cout << "Eig: " << f*eig << endl;
//    cout << "Eig: " << eigf1 << endl;
    cout << "Eig: " << eigf2 << endl;
    cout << "Test integral: " << f_singlet_integral(T) << endl;
    return 0;

    // Testing to confirm Eigen didn't mess up the first vector at least
    //Eigenvector first_vec = vecs.col(0).real();
    //if ( ( P*first_vec 
    //            - vals(0).real() * first_vec ).norm() > 0.00001 ) {
    //    cout << "Matrix Decomposition Failed\n";
    //}
    //double mag = first_vec.transpose() * first_vec;
    //if ( fabs(mag - 1.0) > 0.01 ) cout << "Eigenvector not normalized\n";
    //vals = vals * f;
>>>>>>> origin/main

    cout << "Saving Potential and Susceptibility Functions\n";
    //save_potential_vs_q(FS, P, "potential.dat");
    //if (cube.size() != 0 ) 
    //    save_chi_vs_q(cube, FS, "chi.dat");














/* 
 * ========================================================================================
 * =========================== EIGENVALUE/VECTOR SORTING AND SAVING  ======================
   ========================================================================================
 */
    // Sort solutions with highest eigenvalue/eigenvector pair first
<<<<<<< HEAD
    cout << "Sorting Eigenvectors..." << endl;
    Eigenvector* solutions = lapack_solutions;
    //sort(solutions.rbegin(), solutions.rend());
    vector_to_wave(FS, solutions);
=======
    std::vector<Eigenvector> solutions;
    //solutions = combine_eigs_and_vecs(vals.real(), vecs.real());
    sort(solutions.rbegin(), solutions.rend());
    vector_to_wave(freq_FS[0], solutions);
    //for (int i = 0; i < 6; i++) cout << solutions[0].vec(i) << " "; cout << endl;
    //vector_to_wave(FS, solutions);
    //for (int i = 0; i < 6; i++) cout << solutions[0].vec(i) << " "; cout << endl;
>>>>>>> origin/main
    
    // Defining file name based on cfg (config)
    cout << "Saving Eigenvectors..." << endl;
    std::ostringstream out;
    out.precision(1);
    out << std::fixed << "../data/" + potential_name << dim << "D" 
        << "_mu=" << mu << "_U=" << U << "_wD=" << wc 
        << "_n=" << n << ".dat";
    string file_name = std::move(out).str();

    // Save file in cartesian coordinates for the sake of plotting easier
    save(file_name, T, FS, solutions);
    delete [] solutions;














    return 0;
}
