#include "renormalization.hpp"
#include "src/config/load/cpp_config.hpp"
#include "src/hamiltonian/models/interaction.hpp"
#include "src/objects/vec.hpp"
#include "fields.hpp"
#include <filesystem>

using namespace std;
namespace fs = std::filesystem;

Renormalization::Renormalization() {
    string filename = outdir + prefix + "_renormalization." + filetype;
    if (fs::exists(filename) and automatic_file_read) {
        field = Field_R(filename);
        file_found = true;
    } else {
        file_found = false;
        printf("Renormalization File not found. Defaulting to specified interaction\n");
    }
}

float Renormalization::operator()(Vec k, string label1,
                                  string label2) {
    if (!file_found) return qp_weight;
    return field(k);
}

