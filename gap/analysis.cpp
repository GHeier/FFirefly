/**
 * @file analysis.cpp
 *
 * @brief Contains functions for analyzing the results of the simulation.
 *
 * @author Griffin Heier
 */
#include <iomanip>
#include <iostream>
#include <fstream>

#include <math.h>
#include <vector>
#include <memory>

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
        init_config(mu, U, t, tn, wc, mu, U, t, new_tn, wc);
        vector<Vec> FS = tetrahedron_method(e_base_avg, Vec(0,0,0), mu);
        double T = get_Tc(FS);
        file << tn << " " << T << endl;
    }
}

