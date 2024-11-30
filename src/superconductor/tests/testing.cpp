#include <iostream>
#include <fstream>

#include <cstring>
#include <string>
#include <vector>
#include <complex>

#include <unordered_map>
#include <boost/functional/hash.hpp>
#include <gsl/gsl_integration.h>
#include <omp.h>

#include "../gap/band_structure.h"
#include "../gap/cfg.h"
#include "../gap/eigenvec.hpp"
#include "../gap/frequency_inclusion.hpp"
#include "../gap/matrix.hpp"
#include "../gap/save_data.h"
#include "../gap/susceptibility.h"
#include "../gap/vec.h"
#include "../gap/matrix_creation.h"
#include "../gap/linear_algebra.h"
#include "../gap/solver.h"

#include "integral_tests.h"
#include "matsubara_tests.h"
#include "interpolate_test.h"
#include "plot.h"

using namespace std;

void integral_convergence(float T) {
    int n = 50;
    float mag = 0.04*M_PI;
    Vec q(mag, mag, mag);
    float c2 = integrate_susceptibility(q, T, mu, 0.5, 800);
    for (int i = 0; i < 10; i++) {
        float c1 = integrate_susceptibility(q, T, mu, 0.5, n);
        cout << "Points: " << n << " , Results: " << c1 << " " << c2 << endl;
        n += 50;
    }
}

void mu_to_n() {
    for (int i = 0; i < 5000; i++) {
        float newmu = -0.4 + 0.2/1000.0 * i;
        cout << "Mu: " << newmu << 
            " n: " << 2 * integrate_susceptibility(Vec(0,0,0), 0, newmu, 0, 1000) << endl;
    }
}

void eig_with_freq(float cutoff) {
    ofstream file("eig_freq.dat", std::ios_base::app);
    float T = 0.01;
    init_config(mu, U, t, tn, wc, mu, U, t, tn, cutoff);
    
    vector<vector<Vec>> freq_FS = freq_tetrahedron_method(mu);
    vector<Vec> FS = get_FS(mu);
    printf("FS created\n");

    unordered_map <float, vector<vector<vector<float>>>> cube_freq_map;
    if (potential_name.find("scalapino") != string::npos)
        cube_freq_map = chi_cube_freq(T, mu);
    printf("Cube Created\n");

    int size = 0; for (auto x : freq_FS) size += x.size();
    printf("Size: %i\n", size);

    Matrix P(size); create_P_freq(P, freq_FS, T, cube_freq_map); 
    printf("Matrix Created\n");
    Matrix P2(FS.size()); create_P(P2, FS, T, cube_freq_map);
    printf("Matrix Created\n");

    printf("Sample P values: %f %f %f %f\n", P(0,0), P(0,1), P(1,0), P(1,1));
    printf("Sample P2 values: %f %f %f %f\n", P2(0,0), P2(0,1), P2(1,0), P2(1,1));
    
    float f = f_singlet_integral(T);

    vector<Eigenvector> answers = power_iteration(P, 0.001);
    vector<Eigenvector> answers2 = power_iteration(P2, 0.001);

    float eig = answers[answers.size() - 1].eigenvalue;
    float eig2 = answers2[answers2.size() - 1].eigenvalue;

    cout << "w=0, w>0 Eigs: " << f*eig2 << " " << eig << endl;
    file << cutoff << " " << f*eig2 << " " << eig << endl;
}

void test_cube_map() {
    unordered_map<float, vector<vector<vector<float>>> > cube_freq_map;
    cube_freq_map = chi_cube_freq(0.25, -1.2);
    auto cube = chi_cube(0.25, -1.2, 0.0, "Cube 1/1");
    for (int i = 0; i < cube.size(); i++) {
        for (int j = 0; j < cube[i].size(); j++) {
            for (int k = 0; k < cube[i][j].size(); k++) {
                if (fabs(cube[i][j][k] - cube_freq_map.at(0.0)[i][j][k]) > 0.0001) {
                    printf("Cube values: %f %f\n", cube[i][j][k], cube_freq_map.at(0.0)[i][j][k]);
                }
                if (i == 0) cout << cube[i][j][k] << " " << cube_freq_map.at(0.0)[i][j][k] << endl;
            }
        }
    }
}

int main() {
    
    int num_procs = omp_get_num_procs();
    omp_set_num_threads(num_procs - 1);

    plot_interpolated_chi();
    return 0;
}
