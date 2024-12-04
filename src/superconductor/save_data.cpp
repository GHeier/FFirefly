/**
 * @file save_data.cpp
 *
 * @brief File contains different functions to save calculated data. save() is for Gap Function
 *
 * @author Griffin Heier
 */

#include <iostream>
#include <fstream>
#include <cassert>
#include <string>
#include <math.h>
#include <memory>

#include "../objects/vec.h"
#include "../config/load/cpp_config.h"
#include "../hamiltonian/band_structure.h"
#include "../objects/matrix.hpp"
#include "../objects/eigenvec.hpp"
#include "frequency_inclusion.hpp"

using namespace std;

void save(string file_name, float T, vector<Vec> FS, Eigenvector* solutions) {
    cout << "Saving Eigenvectors..." << endl;
    ofstream file(file_name);
    //file.open(file_name);
    if (dim == 3) file << "x y z ";
    if (dim == 2) file << "x y ";

    for (unsigned int i = 0; i < num_eigenvalues_to_save; i++) 
        file << solutions[i].eigenvalue << " ";
    file << endl;
    for (unsigned int i = 0; i < FS.size(); i++) {
        Vec q = FS[i]; 
        file << q; 
        for (unsigned int j = 0; j < num_eigenvalues_to_save; j++) {
            file << solutions[j].eigenvector[i] << " ";
        }
        file << endl;
    }
    file.close();
    printf("Saved Gap to %s\n", file_name.c_str());
}

void save_with_freq(string file_name, float T, vector<vector<Vec>> &freq_FS, Eigenvector* solutions) {
    printf("Saving Eigenvectors with Frequency...\n");
    ofstream file(file_name);
    if (dim == 3) file << "x y z ";
    if (dim == 2) file << "x y ";

    int size = matrix_size_from_freq_FS(freq_FS);
    for (unsigned int i = 0; i < num_eigenvalues_to_save; i++) {
        file << solutions[i].eigenvalue << " ";
    }
    file << endl;
    int ind = 0;
    printf("Size: %d\n", size);
    printf("Eigenvector size: %d\n", solutions[0].size);
    for (unsigned int i = 0; i < freq_FS.size(); i++) {
        for (unsigned int j = 0; j < freq_FS[i].size(); j++) {
            Vec q = freq_FS[i][j]; 
            file << q; 
            for (unsigned int k = 0; k < num_eigenvalues_to_save; k++) {
                if (ind >= size) {
                    printf("Breaking\n");
                }
                file << solutions[k][ind] << " ";
            }
            ind++;
            file << endl;
        }
    }
    file.close();
    printf("Saved Gap (w-dependent) to %s\n", file_name.c_str());
}

void save_FS(vector<Vec> FS) {
    printf("Saving Fermi Surface...\n");
    ofstream file;
    file.open("FS.dat");
    for (Vec k : FS)
        file << k << endl;
    file.close();
    printf("Saved Fermi Surface to FS.dat\n");
}

void save_potential_vs_q(vector<Vec> &FS, Matrix &P, string filename) {
    printf("Saving Potential...\n");
    ofstream file(filename);
    for (unsigned int i = 0; i < FS.size(); i++) {
        Vec k1 = FS[i];
        for (unsigned int j = 0; j < FS.size(); j++) {
            Vec k2 = FS[j];
            float V = -P(i,j) * pow(vp(k2)/k2.area,0.5) * pow(vp(k1)/k1.area,0.5);
            Vec q = k1 - k2; 
            file << q << V << endl;
        }
    }
    file.close();
    printf("Saved Potential to %s\n", filename.c_str());
}

void save_chi_vs_q(const vector<vector<vector<float>>> &cube, vector<Vec> &FS, string filename) {
    printf("Saving Susceptibility...\n");
    ofstream file(filename);
    for (unsigned int i = 0; i < FS.size(); i++) {
        Vec k1 = FS[i];
        for (unsigned int j = 0; j < FS.size(); j++) {
            Vec k2 = FS[j];
            Vec q = k1 - k2;
            float chi = calculate_chi_from_cube(cube, q);
            file << q << chi << endl;
        }
    }
    file.close();
    printf("Saved Susceptibility to %s\n", filename.c_str());
}

void save_matsubara_cube(const MatCube &cube, float wmin, float wmax, string filename) {
    printf("Saving Matsubara Cube...\n");
    ofstream file(filename);
    file << "Size of each dimension, then min&max values for each dimension\n" << 
        cube.cube.size() << " " 
        << cube.cube[0].size() << " " 
        << cube.cube[0][0].size() << " " 
        << cube.cube[0][0][0].size() << " "
        << 0 << " " << 2*k_max << " "
        << 0 << " " << 2*k_max << " "
        << 0 << " " << 2*k_max << " "
        << wmin << " " << wmax << endl;

    int mx = cube.cube.size();
    int my = cube.cube[0].size();
    int mz = cube.cube[0][0].size();
    int mw = cube.cube[0][0][0].size();

    for (unsigned int i = 0; i < mx; i++) {
        for (unsigned int j = 0; j < my; j++) {
            for (unsigned int k = 0; k < mz; k++) {
                Vec q(2*k_max*i/(mx-1) , 2*k_max*j/(my-1), 2*k_max*k/(mz-1));
                if (mz == 1) q = Vec(2*k_max*i/(mx-1) , 2*k_max*j/(my-1), 0);
                for (unsigned int l = 0; l < mw; l++) {
                    float w = wmin + (wmax - wmin) * l / (mw - 1);
                    file << q; if (dim == 2) file << "0 ";
                    file << w << " (" << cube.cube[i][j][k][l].real() 
                        << ", " << cube.cube[i][j][k][l].imag() << ") \n";
                    if (w == 0 and cube.cube[i][j][k][l].real() == 0) {
                        //printf("Zero at %d %d %d %d\n", i, j, k, l);
                        //printf("q: %f %f %f\n", q.vals[0], q.vals[1], q.vals[2]);
                    }
                }
                file << endl;
            }
        }
    }
    file.close();
    printf("Saved Matsubara Cube to %s\n", filename.c_str());
}
