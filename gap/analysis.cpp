#include <iomanip>
#include <iostream>
#include <fstream>

#include <math.h>
#include <vector>

#include "cfg.h"
#include "fermi_surface.h"
#include "potential.h"
#include "vec.h"
#include "matrix.hpp"
#include "eigenvec.hpp"
#include "save_data.h"
#include "utilities.h"

using namespace std;

void plot_hopping_integral_ratios() {
    ofstream file("crit_temps.dat");
    for (int i = 0; i < 10; i++) {
        double new_tn = (0.4-0.1)*i/10 + 0.1;
        init_config(mu, U, t, tn, w_D, mu, U, t, new_tn, w_D);
        vector<Vec> FS = tetrahedron_method(e_base_avg, Vec(0,0,0), mu);
        double T = get_Tc(FS);
        file << tn << " " << T << endl;
    }
}

void plot_coupling_constants() {
    ofstream file("coupling_constants.dat");
    for (int i = 0; i < 10; i++) {
        double new_tn = (0.4-0.1)*i/10 + 0.1;
        init_config(mu, U, t, tn, w_D, mu, U, t, new_tn, w_D);
        vector<Vec> FS = tetrahedron_method(e_base_avg, Vec(0,0,0), mu);
        double T = get_Tc(FS);
        file << tn << " " << T << endl;
    }
}
