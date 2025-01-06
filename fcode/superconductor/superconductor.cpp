/**
 * Main file for the program
 *
 * This file finds the Fermi Surface(s), calculates the critical temperature, finds the 
 * pairing symmetry, and saves the Gap functions to a file.
 *
 * Author: Griffin Heier
 */
#include <Python.h>

#include <stdio.h>
#include <iomanip>
#include <iostream>
#include <string>

#include <algorithm>
#include <omp.h>
#include <cassert>


#include "../config/load/cpp_config.hpp"
#include "../config/load/py_interface.h"
#include "../config/load/c_config.h"
#include "utilities.hpp"
#include "matrix_creation.hpp"
#include "../algorithms/linear_algebra.hpp"
#include "solver.hpp"
#include "../response/susceptibility.hpp"
#include "save_data.hpp"
#include "../objects/vec.hpp"
#include "../objects/matrix.hpp"
#include "../objects/eigenvec.hpp"
#include "../hamiltonian/band_structure.hpp"
#include "../hamiltonian/interaction.hpp"
#include "cfg.hpp"
#include "superconductor.hpp"

using namespace std;

extern "C" void superconductor_wrapper() {
    printv("Running superconductor_wrapper\n");
    if (method == "bcs")
        bcs();
    else if (method == "eliashberg")
        eliashberg();
    else
        cout << "Method " << method << " not recognized" << endl;
}

void bcs() {
    cout << "Calculating Fermi Surface..." << endl;
    load_cpp_cfg();
 
    vector<vector<Vec>> freq_FS;
    vector<Vec> FS;
    if (not FS_only) {
        freq_FS = freq_tetrahedron_method(mu);
        FS = freq_FS[(l+1)/2 - 1];
    }
    else {
        FS = get_FS(mu);
    }

    cout << "Number of points along Fermi Surface: " << FS.size() << endl;
    float DOS = get_DOS(FS);
    assert(FS.size() > 10);
    save_FS(FS);

    float T = 0.25;
    cout << setprecision(10);
    //cout << coupling_calc(FS, T) << endl;
    //T = 0.065;
    //T = get_Tc(FS);
    printf("Temperature: %.5f \n", T);

    // Calculates the susceptibility matrix if it's going to be used in the potential
    // Otherwise it's passed as empty

    load_chi("chi_mesh_dynamic.dat");

    int m_size = FS.size();
    if (not FS_only) m_size = matrix_size_from_freq_FS(freq_FS);

    Matrix P(m_size);
    create_P(P, FS);
    float f = f_singlet_integral(T);
    cout << "F-integral value: " << f << endl;

    cout << "Finding Eigenspace..." << endl;
    Eigenvector *solutions = new Eigenvector[num_eigenvalues_to_save];
    lapack_hermitian_diagonalization(P, solutions);

    // Sort solutions with highest eigenvalue/eigenvector pair first
    cout << "Sorting Eigenvectors..." << endl;
    sort(solutions, solutions + num_eigenvalues_to_save, descending_eigenvalues);
    cout << "Sorted Eigenvectors\n";
    if(FS_only) vector_to_wave(FS, solutions);
    else freq_vector_to_wave(freq_FS, solutions);
    
    // Defining file name based on config/load/cpp_config (config)
    cout << "Saving Eigenvectors..." << endl;
    string file_name = get_SC_filename();
    cout << "File Name: " << file_name << endl;

    // Save file in cartesian coordinates for the sake of plotting easier
    if (FS_only) save(file_name, T, FS, solutions);
    else save_with_freq(file_name, T, freq_FS, solutions);
    cout << "Eigenvectors Saved\n";
    delete [] solutions;
}

void eliashberg() {
    string folder = "superconductor/";
    string filename = "eliashberg";
    string function = "eliashberg";
    call_python_func(folder.c_str(), filename.c_str(), function.c_str());
}
