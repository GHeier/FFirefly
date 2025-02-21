#include <iostream>
#include <math.h>
#include <vector>
#include <string>

#include <omp.h>
#include <boost/functional/hash.hpp>
#include <unordered_map>

#include "interaction.hpp"
#include "../hamiltonian/band_structure.hpp"
#include "../objects/vec.hpp"
#include "../config/load/cpp_config.hpp"
//#include "../response/susceptibility.hpp"

using namespace std;

// Potential functions
__attribute__((visibility("default")))
float V(const Vec q, float w, string spin1, string spin2) {
    if (interaction == "const") 
        return potential_const();
    if (interaction == "phonon_coulomb") 
        return phonon_coulomb(q);
    if (interaction == "test") 
        return potential_test(q, spin1, spin2);
    else {
        cout << "Unknown Potential Function: " << interaction << endl;
        exit(1);
    }
}


float potential_const() {
    return -1;
}

float potential_test(Vec k1, string spin1, string spin2) {
    Vec q1 = k1;
    return -1.0;
}

float phonon_coulomb(Vec q) {
    float qx = q(0);
    float Vp = 1.0/3.0;
    if (q.norm() != 0) {
        Vp = 1/(1+2*qx*qx / pow(q.norm(), 2));
    }
    float Vc = 1 / (1 + q.norm());
    return Vp + Vc;
}

/* ======================================================================
* ======================== C Interaction Functions ========================
*/

float Vs_c(double k1_c[3], double k2_c[3], const char* spin1_c, const char* spin2_c) {
    Vec k1 = Vec(k1_c[0], k1_c[1], k1_c[2]);
    Vec k2 = Vec(k2_c[0], k2_c[1], k2_c[2]);
    string spin1(spin1_c);
    string spin2(spin2_c);
    if (interaction == "const") 
        return potential_const();
    if (interaction == "test") 
        return potential_test(k1-k2, spin1, spin2);
    else {
        cout << "Unknown Potential Function: " << interaction << endl;
        exit(1);
    }
}


float V_c(double k1_c[3], double k2_c[3]) {
    Vec k1 = Vec(k1_c[0], k1_c[1], k1_c[2]);
    Vec k2 = Vec(k2_c[0], k2_c[1], k2_c[2]);
    string spin1 = "up";
    string spin2 = "up";
    if (interaction == "const") 
        return potential_const();
    if (interaction == "test") 
        return potential_test(k1-k2, spin1, spin2);
    else {
        cout << "Unknown Potential Function: " << interaction << endl;
        exit(1);
    }
}
