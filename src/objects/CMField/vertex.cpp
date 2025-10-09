#include "vertex.hpp"
#include "../../config/load/cpp_config.hpp"
#include "../../hamiltonian/interaction.hpp"
#include "../vec.hpp"
#include "fields.hpp"
#include <filesystem>

#include <complex>

using namespace std;
namespace fs = std::filesystem;

Vertex::Vertex() {
    string filename = outdir + prefix + "_vertex." + filetype;
    if (fs::exists(filename) and automatic_file_read) {
        field = Field_C(filename);
        file_found = true;
    } else {
        file_found = false;
        printf("Vertex File not found. Defaulting to specified interaction\n");
    }
}

complex<float> Vertex::operator()(Vec k, float w, string label1,
                                  string label2) {
    if (!file_found)
        return complex<float>(V(k, imag(w), label1, label2), 0);
    return field(k, w);
}

complex<float> Vertex::operator()(Vec k, complex<float> w, string label1,
                                  string label2) {
    return operator()(k, imag(w), label1, label2);
}
