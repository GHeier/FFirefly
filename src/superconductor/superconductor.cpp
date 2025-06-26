/**
 * Main file for the program
 *
 * This file finds the Fermi Surface(s), calculates the critical temperature,
 * finds the pairing symmetry, and saves the Gap functions to a file.
 *
 * Author: Griffin Heier
 */
#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <string>

#include <algorithm>
#include <cassert>
#include <omp.h>

#include "../algorithms/linear_algebra.hpp"
#include "../config/load/cpp_config.hpp"
#include "../config/load/jl_interface.h"
#include "../hamiltonian/band_structure.hpp"
#include "../objects/eigenvec.hpp"
#include "../objects/matrix.hpp"
#include "../objects/vec.hpp"
#include "cfg.hpp"
#include "matrix_creation.hpp"
#include "save_data.hpp"
#include "solver.hpp"
#include "superconductor.hpp"
#include "utilities.hpp"

using namespace std;

void bcs() {
    cout << "Calculating Fermi Surface..." << endl;
    load_cpp_cfg();

    vector<vector<Vec>> freq_FS;
    vector<Vec> FS;
    if (not FS_only) {
        freq_FS = freq_tetrahedron_method(mu);
        FS = freq_FS[(l + 1) / 2 - 1];
    } else {
        FS = get_FS(mu);
    }

    cout << "Number of points along Fermi Surface: " << FS.size() << endl;
    //float DOS = get_DOS(FS);
    //printf("DOS: %f\n", DOS);
    assert(FS.size() > 10);
    //save_FS(FS);

    float T = Temperature;
    cout << setprecision(10);
    // cout << coupling_calc(FS, T) << endl;
    // T = 0.065;
    // T = get_Tc(FS);
    printf("Temperature: %.5f \n", T);

    float renorm = 1.0;
    if (FS_only)
        renorm = get_renormalization(FS);
    else
        renorm = get_renormalization_off_FS(freq_FS);
    printf("renorm = %f\n", renorm);

    // Calculates the susceptibility matrix if it's going to be used in the
    // potential Otherwise it's passed as empty

    int m_size = FS.size();
    if (not FS_only)
        m_size = matrix_size_from_freq_FS(freq_FS);

    Matrix P(m_size);
    if (FS_only)
        create_P(P, FS);
    else {
        create_P_freq(P, freq_FS, T);
    }
    float f = f_singlet_integral(T);
    cout << "F-integral value: " << f << endl;

    vector<float> proj_eigs = matrix_projections(FS, P, renorm);

    //Eigenvector initial_guess(P.size, true);
    //for (int i = 0; i < P.size; i++) {
    //    Vec k = FS[i];
    //    float proj = cos(k(0)) - cos(k(1));
    //    initial_guess.eigenvector[i] = proj;
    //}
    Eigenvector top_gap = power_iteration(P);
    //Eigenvector* temp = new Eigenvector[1];
    //temp[0] = top_gap;
    printf("Max Power Iteration eigenvalue: %f\n", top_gap.eigenvalue / (1 + renorm));
    //vector<Eigenvector> top_gaps = power_iteration(P, 0.01);
    //for (Eigenvector x : top_gaps) {
    //    cout << "eig: " << x.eigenvalue << endl;
    //}

    cout << "Finding Eigenspace..." << endl;
    Eigenvector *solutions = new Eigenvector[num_eigenvalues_to_save];
    lapack_hermitian_diagonalization(P, solutions);

    // Sort solutions with highest eigenvalue/eigenvector pair first
    cout << "Sorting Eigenvectors..." << endl;
    sort(solutions, solutions + num_eigenvalues_to_save,
         descending_eigenvalues);

    printf("Max Diagonalized eigenvalue: %f\n", solutions[0].eigenvalue / (1 + renorm));

    if (FS_only) {
        double calc_Tc = get_Tc_FS_only(solutions[0].eigenvalue);
        printf("Calculated Tc: %.5f\n", calc_Tc);
    }
    cout << "Sorted Eigenvectors\n";
    if (FS_only)
        vector_to_wave(FS, solutions);
    else
        freq_vector_to_wave(freq_FS, solutions);

    auto d_x2_y2 = [](Vec k) { return cos(k(0)) - cos(k(1)); };
    float norm = 0;
    for (int i = 0; i < FS.size(); i++) {
        Vec k = FS[i];
        norm += d_x2_y2(k) * d_x2_y2(k);
    }
    for (int i = 0; i < FS.size(); i++) {
        Vec k = FS[i];
        //cout << d_x2_y2(k) / pow(norm, 0.5) << ", " << solutions[0].eigenvector[i] << endl;
    }

    // Defining file name based on config/load/cpp_config (config)
    cout << "Saving Eigenvectors..." << endl;
    string file_name = get_SC_filename();
    file_name = outdir + prefix + "_gap." + filetype;

    // Save file in cartesian coordinates for the sake of plotting easier
    if (FS_only)
        save(file_name, T, FS, solutions);
        //save(file_name, T, FS, temp);
    else
        save_with_freq(file_name, T, freq_FS, solutions);
    cout << "Eigenvectors Saved\n";
    delete[] solutions;
}

void eliashberg() {
    string folder = "superconductor/";
    string filename = "eliashberg";
    string module = "Eliashberg";
    string function = "eliashberg_node";
    // call_python_func(folder.c_str(), filename.c_str(), function.c_str());
    call_julia_func(folder.c_str(), filename.c_str(), module.c_str(),
                    function.c_str());
}

void linearized_eliashberg() {
    string folder = "superconductor/";
    string filename = "linearized_eliashberg";
    string module = "Linearized_Eliashberg";
    string function = "eigenvalue_computation";
    // call_python_func(folder.c_str(), filename.c_str(), function.c_str());
    call_julia_func(folder.c_str(), filename.c_str(), module.c_str(),
                    function.c_str());
}

void debug() {
    string folder = "superconductor/";
    string filename = "lin_eliashberg_surface";
    string module = "Linearized_Eliashberg_Surface";
    string function = "eigenvalue_computation";
    // call_python_func(folder.c_str(), filename.c_str(), function.c_str());
    call_julia_func(folder.c_str(), filename.c_str(), module.c_str(),
                    function.c_str());
}
