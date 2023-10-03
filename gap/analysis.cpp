#include <iomanip>
#include <iostream>
#include <fstream>

#include <math.h>
#include <vector>

#include "cfg.h"
#include "fermi_surface.h"
#include "potential.h"
#include "vec.h"
#include "save_data.h"
#include "utilities.h"

using namespace std;

void plot_hopping_integral_ratios() {
    ofstream file("crit_temps.dat");
    for (int i = 0; i < 10; i++) {
        double new_tn = (0.4-0.1)*i/10 + 0.1;
        init_config(mu, U, t, tn, mu, U, t, new_tn);
        vector<Vec> FS = tetrahedron_method(mu);
        double T = get_Tc(FS);
        file << tn << " " << T << endl;
    }
}
