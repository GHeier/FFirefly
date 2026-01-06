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

#include "src/algorithms/linear_algebra.hpp"
#include "src/config/load/cpp_config.hpp"
#include "src/config/load/py_interface.h"
#include "src/config/load/jl_interface.h"
#include "src/hamiltonian/models/band_structure.hpp"
#include "src/objects/CMField/fields.hpp"
#include "src/objects/eigenvec.hpp"
#include "src/objects/matrix.hpp"
#include "src/objects/vec.hpp"
#include "cfg.hpp"
#include "matrix_creation.hpp"
#include "save_data.hpp"
#include "solver.hpp"
#include "superconductor.hpp"
#include "utilities.hpp"

using namespace std;

float bcs() {
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

    float renorm = 0.0;
    if (FS_only)
        renorm = get_renormalization(FS);
    else
        renorm = get_renormalization_off_FS(freq_FS);
    double onesum = 0.0;
    Field_C der_sigma(outdir + prefix + "_renormalization.h5");
    for (Vec x : FS) {
        onesum += (x.area / vp(x.n, x) * real(der_sigma(x, 0.0)));
    }
    printf("Average analytic dSigma/dw on FS: %f\n", onesum / (pow(2 * M_PI, dim)));
    printf("lambda_z = %f\n", renorm);

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

    printf("\n");
    vector<float> proj_eigs = matrix_projections(FS, P, renorm);
    printf("\n");

    //Eigenvector initial_guess(P.size, true);
    //for (int i = 0; i < P.size; i++) {
    //    Vec k = FS[i];
    //    float proj = cos(k(0)) - cos(k(1));
    //    initial_guess.eigenvector[i] = proj;
    //}
    Eigenvector *solutions = new Eigenvector[num_eigenvalues_to_save];
    //if (method == "power_iteration") {
        printf("Performing Power Iteration\n");
        Eigenvector top_gap = power_iteration(P);
        printf("Max Power Iteration eigenvalue: %f\n", top_gap.eigenvalue);
        printf("Max Effective eigenvalue: %f\n", top_gap.eigenvalue / (1 + renorm));
        printf("Max Power Iteration Eigenvalue with T included: %f\n", top_gap.eigenvalue / (1 + renorm) * f);
        printf("Max Effective Eigenvalue with T included: %f\n", top_gap.eigenvalue / (1 + renorm) * f);
        solutions[0] = top_gap;
    //}
    //else if (method == "diagonalization") {
    //    cout << "Finding Eigenspace..." << endl;
    //    vector<Eigenvector> temp_solutions = lapack_hermitian_diagonalization(P);
    //    for (int i = 0; i < num_eigenvalues_to_save; i++) {
    //        solutions[i] = temp_solutions[i];
    //    }
    //    temp_solutions.clear();

    //    // Sort solutions with highest eigenvalue/eigenvector pair first
    //    cout << "Sorting Eigenvectors..." << endl;
    //    sort(solutions, solutions + num_eigenvalues_to_save,
    //        descending_eigenvalues);

    //    printf("Max Diagonalized eigenvalue: %f\n", solutions[0].eigenvalue / (1 + renorm));
    //    printf("Eigenvalue with T included: %f\n", solutions[0].eigenvalue / (1 + renorm) * f);
    //}
    //else {
    //    printf("No method provided. Exiting\n");
    //    exit(0);
    //}
    printf("\n");
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
        //save(file_name, T, FS, temp);
    else
        save_with_freq(file_name, T, freq_FS, solutions);
    cout << "Eigenvectors Saved\n";
    delete[] solutions;
    return top_gap.eigenvalue;
}

