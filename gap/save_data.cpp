/**
 * @file save_data.cpp
 *
 * @brief File contains different functions to save calculated data. save() is for Gap Function
 *
 * @author Griffin Heier
 */

#include <iomanip>
#include <sstream>
#include <filesystem>
#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>
#include <string>
#include <math.h>
#include <unordered_map>
#include <memory>

#include "calculations.h"
#include "frequency_inclusion.hpp"
#include "potential.h"
#include "susceptibility.h"
#include "vec.h"
#include "cfg.h"
#include "matrix.hpp"
#include "eigenvec.hpp"

using std::endl;
using std::cout;
using std::vector;
using std::string;

namespace fs = std::filesystem;

void save(string file_name, float T, vector<Vec> FS, Eigenvector* solutions) {
    std::ofstream file(file_name);
    //file.open(file_name);
    if (dim == 3) file << "x y z ";
    if (dim == 2) file << "x y ";

    for (unsigned int i = 0; i < num_eigenvalues_to_save; i++) 
        file << solutions[i].eigenvalue << " ";
    file << endl;
    for (unsigned int i = 0; i < FS.size(); i++) {
        Vec q = FS[i]; 
        if(not q.cartesian) 
            q.to_cartesian();
        file << q; 
        for (unsigned int j = 0; j < num_eigenvalues_to_save; j++) {
            file << solutions[j].eigenvector[i] << " ";
        }
        file << endl;
    }
}

void save_with_freq(string file_name, float T, vector<vector<Vec>> &freq_FS, Eigenvector* solutions) {
    std::ofstream file(file_name);
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
            if(not q.cartesian) 
                q.to_cartesian();
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
}





void save_FS(const vector<Vec>& FS, int iteration, float mu) {
    // completely unnecessary check for existence of FS directory
    if (!fs::exists("FS")) {
        fs::create_directory("FS");
    }

    // cleaning out fs directory on the first iteration
    if (iteration == 0) {
        for (const auto& entry : fs::directory_iterator("FS")) {
            try {
                fs::remove(entry);
            } catch (const fs::filesystem_error& e) {
                cerr << "Error deleting file: " << entry.path() << " - " << e.what() << endl;
            }
        }
    }

    // naming file. include mu for graphing func later

    int mu_prec = 5; // float precision for mu n filename

    ostringstream filename;
    filename << "FS/FS_prec_" << mu_prec << "_mu_" << fixed << setprecision(mu_prec) << mu << "_iter_" << iteration << ".dat";


    std::ofstream file(filename.str());
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename.str() << std::endl;
        return;
    }

    // write data
    for (Vec k : FS) {
        file << k << endl;
    }
    
    file.close();
}


void save_FS(vector<Vec> FS) {
    // completely unnecessary check for existence of FS directory
    if (!fs::exists("FS")) {
        fs::create_directory("FS");
    }

    // cleaning out fs directory on the first iteration
    for (const auto& entry : fs::directory_iterator("FS")) {
            try {
                fs::remove(entry);
            } catch (const fs::filesystem_error& e) {
                cerr << "Error deleting file: " << entry.path() << " - " << e.what() << endl;
            }
        }

    std::ofstream file;
    file.open("FS/FS.dat"); // for now manually delete the FS files on succesive iterations
    
    // writing file
    for (Vec k : FS)
        file << k << endl;
}

void save_potential_vs_q(vector<Vec> &FS, Matrix &P, string filename) {
    ofstream file(filename);
    for (unsigned int i = 0; i < FS.size(); i++) {
        Vec k1 = FS[i];
        for (unsigned int j = 0; j < FS.size(); j++) {
            Vec k2 = FS[j];
            float V = -P(i,j) * pow(vp(k2)/k2.area,0.5) * pow(vp(k1)/k1.area,0.5);
            Vec q = k1 - k2; if (q.cartesian == false) q.to_cartesian();
            file << q << V << endl;
        }
    }
}

void save_chi_vs_q(const vector<vector<vector<float>>> &cube, vector<Vec> &FS, string filename) {
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
}
