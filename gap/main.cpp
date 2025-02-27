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
#include "save_data_recursive.h"
#include "vec.h"
#include "matrix.hpp"
#include "eigenvec.hpp"
#include "utilities.h"

using std::string;
float get_max_eigenvalue() {
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
    vector<Vec> layer = tetrahedron_method(e_base_avg, Vec(0,0,0), 0);
    
    printf("\n");
    for (int i = 0; i < 5; i++)
    {
	Vec k = layer[i]; // was once vec i believe this was an error
	printf("v_k = %lf\n", vp(k));
    }
    printf("\n");

    save_vp_dir(mu, layer);
    save_V_dir(mu, layer);

    // this is a commented out section but un commenting it out did not fix anything
    /*
    unordered_map <float, vector<vector<vector<float>>>> nothing;
    ofstream vp_file("vp_file.dat");
    ofstream V_file("V_file.dat");
    for (int i = 0; i < layer.size(); i++) {
	    Vec k = layer[i];
	    vp_file << vp(k) << endl;
    }

    for (int i = 0; i < layer.size(); i++) {
	    for (int j = 0; j < layer.size(); j++) {
		    Vec k1 = layer[i];
		    Vec k2 = layer[j];
		    float V_val = V(k1, k2, 0, 0, nothing);
		    V_file << V_val << endl;
	    }
    }
    */
    


    cout << "Number of points along Fermi Surface: " << FS.size() << endl;
    save_FS(layer, iteration, mu);

    //PUT VELOCITY STUFF HERE (I lied)



    float DOS = get_DOS(FS);
    cout << "Density of States: " << DOS << endl;
    //cout << "Ending code early\n";
    //return 0;












/* 
 * ========================================================================================
 * =========================== CRITICAL TEMPERATURE CALCULATION  ==========================
   ========================================================================================
 */
    float T = 0.25000;
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

    // This calculates total size of the matrix
    int size = matrix_size_from_freq_FS(freq_FS);

    //Matrix Pf2(size); 
    //create_P_freq(Pf2, freq_FS, T, cube_freq_map);
    Matrix P(FS.size());
    create_P(P, FS, T, cube_freq_map);

    float f = f_singlet_integral(T);
    cout << "F-integral value: " << f << endl;

    cout << "Finding Eigenspace..." << endl;
    Eigenvector *lapack_solutions = new Eigenvector[num_eigenvalues_to_save];
    lapack_hermitian_diagonalization(P, lapack_solutions);

    cout << "Saving Potential and Susceptibility Functions\n";
    //save_potential_vs_q(FS, P, "potential.dat");
    //if (cube.size() != 0 ) 
    //    save_chi_vs_q(cube, FS, "chi.dat");




/*
 * ========================================================================================
 * ================= FERMI VELOCITY AND INTERACTION POTENTIAL SAVING ======================
   ========================================================================================
 */
    //fermi_velocity_average(mu);
    //interaction_potential_average(mu);








/* 
 * ========================================================================================
 * =========================== EIGENVALUE/VECTOR SORTING AND SAVING  ======================
   ========================================================================================
 */
    // Sort solutions with highest eigenvalue/eigenvector pair first
    cout << "Sorting Eigenvectors..." << endl;
    Eigenvector* solutions = lapack_solutions;
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

    // Save file in cartesian coordinates for the sake of plotting easier
    save(file_name, T, FS, solutions);
    //delete [] solutions;













	return solutions[0].eigenvalue;
    
    // taking averages	
    fermi_velocity_average(mu);
    interaction_potential_average(mu);

}

int iteration;

int main() {

	std::cout << "\nYou are about to run the gap code. Do you want to clear and archive the data files in the fermi_velocity_data and interaction_potential_data directories? (y/n)\n";

	char response;
	
	std::cin >> response;

	response = std::tolower(response);

	if(response == 'y')
	{

		archive_and_clear_data("fermi_velocity_data"); archive_and_clear_data("interaction_potential_data"); // empty and archive old data files
	}

	else
	{
		std::cout << "\nProceeding without deleting old data.\n";
	}

	int num_iters = 20;
	for (int i = 0; i <= num_iters; i++) {
		
		printf("\n\n\n==============================\nBeginning code for iteration %i\n==============================\n", i);

		int iteration = i; 
		
		float new_mu = 0 - 4.0 * i / num_iters;
		
		//mu = new_mu;


		change_global_constant(mu, new_mu);
		float eig = get_max_eigenvalue();
		printf("Iteration %i | Eigenvalue: %.4f\n", i, eig);
		printf("\n\n");
		printf("\n==============================\nEnding  code  for  iteration %i\n==============================\n\n\n\n\n", i);
	}
   	return 0;
}
