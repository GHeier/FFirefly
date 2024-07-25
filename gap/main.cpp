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
    printf("Density of States: %.5f\n", DOS);











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
    unordered_map <float, vector<vector<vector<float>>>> cube_freq_map;
    // Calculates the susceptibility matrix if it's going to be used in the potential
    // Otherwise it's passed as empty
    if (potential_name.find("scalapino") != string::npos) 
        cube_freq_map = chi_cube_freq(T, mu);


    //Matrix Pf2(size); 
    //create_P_freq(Pf2, freq_FS, T, cube_freq_map);
    Matrix P(FS.size());
    create_P(P, FS, T, cube_freq_map);
    FILE *debug_file = fopen("debug_file.txt", "w");
    for (int i = 0; i < FS.size(); i++) {
        for (int j = 0; j < FS.size(); j++) {
            fprintf(debug_file, "P(%d, %d) = %f\n", i, j, P(i,j));
        }
    }

    float f = f_singlet_integral(T);
    cout << "F-integral value: " << f << endl;

    cout << "Finding Eigenspace..." << endl;
    Eigenvector solutions[P.size];
    lapack_diagonalization(P, solutions);

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
    cout << "Sorting Eigenvectors..." << endl;
    //sort(solutions.rbegin(), solutions.rend());
    vector_to_wave(FS, solutions);
    
    // Defining file name based on cfg (config)
    cout << "Saving Eigenvectors..." << endl;
    std::ostringstream out;
    out.precision(1);
    out << std::fixed << "../data/" + potential_name << dim << "D" 
        << "_mu=" << mu << "_U=" << U << "_wD=" << wc 
        << "_n=" << n << ".dat";
    string file_name = std::move(out).str();
    cout << "File Name: " << file_name << endl;

    // Save file in cartesian coordinates for the sake of plotting easier
    save_with_freq(file_name, T, freq_FS, solutions);
    cout << "Eigenvectors Saved\n";
    //delete [] solutions;














    return 0;
}
