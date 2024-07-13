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
    vector<Vec> FS = tetrahedron_method(e_base_avg, Vec(0,0,0), mu);
    if (FS != freq_FS[(l+1)/2 - 1]) {
        cout << "Fermi Surface Calculations Failed\n";
    }
    FS = freq_FS[(l+1)/2 - 1];


    cout << "Number of points along Fermi Surface: " << freq_FS[(l+1)/2 - 1].size() << endl;
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
    printf("Temperature: %.5f \n", T);

/* 
 * ========================================================================================
 * ========================== MATRIX CREATION AND DIAGONALIZATION  ========================
   ========================================================================================
 */
    unordered_map <double, vector<vector<vector<double>>>> cube_freq_map;
    if (potential_name.find("scalapino") != string::npos) 
        cube_freq_map = chi_cube_freq(T, mu);

    int size = 0;
    for (int i = 0; i < freq_FS.size(); i++) {
        size += freq_FS[i].size();
    }
    //Matrix Pf2(size); 
    //create_P_freq(Pf2, freq_FS, T, cube_freq_map);
    Matrix P(FS.size());
    create_P(P, FS, T, cube_freq_map);
    double f = f_singlet_integral(T);
    cout << "F-integral value: " << f << endl;

    cout << "Finding Eigenspace..." << endl;
    vector<Eigenvector> lapack_solutions = lapack_diagonalization(P);
    //vector<Eigenvector> answers = power_iteration(P, 0.001);
    //vector<Eigenvector> answersf2 = power_iteration(Pf2, 0.001);
    vector<Eigenvector> answers = lapack_solutions;
    Eigenvector v1 = answers[answers.size() - 1];
    //Eigenvector v2 = answersf2[answersf2.size() - 1];
    cout << "Eig: " << f*v1.eigenvalue << endl;
    //cout << "Eig: " << v2.eigenvalue << endl;

    // Testing to confirm Eigen didn't mess up the first vector at least
    if ( ( P*v1 - v1 * v1.eigenvalue).norm() > 0.00001 ) {
        cout << "First Matrix Decomposition Failed\n";
    }
    //if ( ( Pf2*v2 - v2 * v2.eigenvalue).norm() > 0.00001 ) {
    //    cout << "Second Matrix Decomposition Failed\n";
    //}
    if ( fabs(v1.norm() - 1.0) > 0.01 ) cout << "Eigenvector not normalized\n";

    cout << "Saving Potential and Susceptibility Functions\n";
    save_potential_vs_q(FS, P, "potential.dat");
    //if (cube.size() != 0 ) 
    //    save_chi_vs_q(cube, FS, "chi.dat");

/* 
 * ========================================================================================
 * =========================== EIGENVALUE/VECTOR SORTING AND SAVING  ======================
   ========================================================================================
 */
    // Sort solutions with highest eigenvalue/eigenvector pair first
    cout << "Sorting Eigenvectors..." << endl;
    std::vector<Eigenvector> solutions;
    solutions.push_back(v1);
    //sort(solutions.rbegin(), solutions.rend());
    vector_to_wave(FS, solutions);
    
    // Defining file name based on cfg (config)
    cout << "Saving Eigenvectors..." << endl;
    std::ostringstream out;
    out.precision(1);
    out << std::fixed << "../data/" + potential_name << dim << "D" 
        << "_mu=" << mu << "_U=" << U << "_wD=" << w_D 
        << "_n=" << n << ".dat";
    string file_name = std::move(out).str();

    // Save file in cartesian coordinates for the sake of plotting easier
    save(file_name, T, FS, solutions);


    return 0;
}
