#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <Eigen/Dense>

#include "calculations.h"
#include "potential.h"
#include "vec.h"
#include "cfg.h"

using namespace Eigen;
using std::endl;
using std::cout;
using std::vector;
using std::string;

void save(string file_name, double T, vector<Vec> k, std::vector<EigAndVec> solutions) {
    std::ofstream file(file_name);
    //file.open(file_name);
    if (dim == 3) file << "x y z ";
    if (dim == 2) file << "x y ";

    for (unsigned int i = 0; i < solutions.size(); i++) file << solutions[i].eig << " ";
    file << endl;
    for (unsigned int i = 0; i < k.size(); i++) {
        Vec q = k[i]; if(not q.cartesian) q.to_cartesian();
        file << q; 
        for (unsigned int j = 0; j < solutions[0].vec.size(); j++) {
            file << solutions[j].vec(i) << " ";
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

void save_potential_vs_q(vector<Vec> &FS, MatrixXd &P, string filename) {
    ofstream file(filename);
    for (unsigned int i = 0; i < FS.size(); i++) {
        Vec k1 = FS[i];
        for (unsigned int j = 0; j < FS.size(); j++) {
            Vec k2 = FS[j];
            double V = -P(i,j) * pow(vp(k2)/k2.area,0.5) * pow(vp(k1)/k1.area,0.5);
            Vec q = k1 - k2; if (q.cartesian == false) q.to_cartesian();
            file << q << V << endl;
        }
    }
}

void save_chi_vs_q(const vector<vector<vector<double>>> &cube, vector<Vec> &FS, string filename) {
    ofstream file(filename);
    for (unsigned int i = 0; i < FS.size(); i++) {
        Vec k1 = FS[i];
        for (unsigned int j = 0; j < FS.size(); j++) {
            Vec k2 = FS[j];
            Vec q = k1 - k2;
            double chi = calculate_chi_from_cube(cube, q);
            file << q << chi << endl;
        }
    }
}
