#include <fstream>
#include "../config/load/cpp_config.hpp"
#include "band_structure.hpp"

extern "C" void save_FS() {
    float E = fermi_energy;
    vector<Vec> FS = get_FS(E);
    printv("Saving Fermi Surface for E = %.2f\n", E);
    ofstream file(outdir + prefix + "_FS.dat");
    if (dimension == 3)
        file << "# x y z n\n";
    else
        file << "# x y n\n";
    for (auto k : FS) {
        file << k << " " << k.n << "\n";
    }
    file.close();
}


