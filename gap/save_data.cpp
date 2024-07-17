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
#include <vector>
#include <string>
#include <math.h>
#include <unordered_map>

#include "calculations.h"
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

<<<<<<< HEAD
void save(string file_name, float T, vector<Vec> FS, Eigenvector* solutions) {
=======
void save(string file_name, double T, vector<Vec> k, std::vector<Eigenvector> solutions) {
>>>>>>> origin/main
    std::ofstream file(file_name);
    //file.open(file_name);
    if (dim == 3) file << "x y z ";
    if (dim == 2) file << "x y ";

<<<<<<< HEAD
    for (unsigned int i = 0; i < FS.size(); i++) 
        file << solutions[i].eigenvalue << " ";
=======
    for (unsigned int i = 0; i < solutions.size(); i++) file << solutions[i].eigenvalue << " ";
>>>>>>> origin/main
    file << endl;
    for (unsigned int i = 0; i < FS.size(); i++) {
        Vec q = FS[i]; 
        if(not q.cartesian) 
            q.to_cartesian();
        file << q; 
<<<<<<< HEAD
        for (unsigned int j = 0; j < FS.size(); j++) {
            file << solutions[j].eigenvector[i] << " ";
=======
        for (unsigned int j = 0; j < solutions[0].size; j++) {
            file << solutions[j][i] << " ";
>>>>>>> origin/main
        }
        file << endl;
    }
}

void save_FS(vector<Vec> FS) {
    std::ofstream file;
    file.open("FS.dat");
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
