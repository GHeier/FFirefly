#include "self_energy.hpp"
#include "../../config/load/cpp_config.hpp"
#include "../../hamiltonian/interaction.hpp"
#include "../vec.hpp"
#include "cmfield.hpp"
#include "fields.hpp"
#include <filesystem>

#include <complex>

using namespace std;
namespace fs = std::filesystem;

Self_Energy::Self_Energy() {
    string filename = outdir + prefix + "_self_energy." + filetype;
    if (fs::exists(filename) and automatic_file_read) {
        field = load_CMField(filename);
        file_found = true;
    } else {
        file_found = false;
        printf("self_energy File not found. Defaulting to specified interaction\n");
    }
}

complex<float> Self_Energy::operator()(Vec k, float w, string label1,
                                  string label2) {
    if (!file_found)
        return complex<float>(0, 0);
    complex<Vec> val = field(k, w);
    return complex<float>(val.real().x, val.imag().x);
}

complex<float> Self_Energy::operator()(Vec k, complex<float> w, string label1,
                                  string label2) {
    return operator()(k, imag(w), label1, label2);
}

