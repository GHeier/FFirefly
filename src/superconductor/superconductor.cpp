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
    float DOS = get_DOS(FS);
    printf("DOS: %f\n", DOS);
    assert(FS.size() > 10);
    save_FS(FS);

    float T = Temperature;
    cout << setprecision(10);
    // cout << coupling_calc(FS, T) << endl;
    // T = 0.065;
    // T = get_Tc(FS);
    printf("Temperature: %.5f \n", T);

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

    cout << "Finding Eigenspace..." << endl;
    Eigenvector *solutions = new Eigenvector[num_eigenvalues_to_save];
    lapack_hermitian_diagonalization(P, solutions);

    // Sort solutions with highest eigenvalue/eigenvector pair first
    cout << "Sorting Eigenvectors..." << endl;
    sort(solutions, solutions + num_eigenvalues_to_save,
         descending_eigenvalues);

    printf("Max eigenvalue: %f\n", solutions[0].eigenvalue);
    if (FS_only) {
        double calc_Tc = get_Tc_FS_only(solutions[0].eigenvalue);
        printf("Calculated Tc: %.5f\n", calc_Tc);
    }
    cout << "Sorted Eigenvectors\n";
    if (FS_only)
        vector_to_wave(FS, solutions);
    else
        freq_vector_to_wave(freq_FS, solutions);

    // Defining file name based on config/load/cpp_config (config)
    cout << "Saving Eigenvectors..." << endl;
    string file_name = get_SC_filename();
    file_name = outdir + prefix + "_gap.dat";

    // Save file in cartesian coordinates for the sake of plotting easier
    if (FS_only)
        save(file_name, T, FS, solutions);
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

void debug() {
    string folder = "superconductor/";
    string filename = "debug";
    string module = "Test";
    string function = "test";
    // call_python_func(folder.c_str(), filename.c_str(), function.c_str());
    call_julia_func(folder.c_str(), filename.c_str(), module.c_str(),
                    function.c_str());
}
