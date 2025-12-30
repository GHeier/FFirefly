#include "self_energy.hpp"
#include "src/config/load/cpp_config.hpp"
#include "src/hamiltonian/models/interaction.hpp"
#include "src/objects/vec.hpp"
#include "fields.hpp"
#include <filesystem>

#include <complex>

using namespace std;
namespace fs = std::filesystem;

Self_Energy::Self_Energy() {
    string filename = outdir + prefix + "_self_energy." + filetype;
    if (fs::exists(filename) and automatic_file_read) {
        field = Field_C(filename);
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
    return field(k, w);
}

complex<float> Self_Energy::operator()(Vec k, complex<float> w, string label1,
                                  string label2) {
    return operator()(k, imag(w), label1, label2);
}

