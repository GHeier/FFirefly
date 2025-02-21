#include "vertex.hpp"
#include "fields.hpp"
#include "../../hamiltonian/interaction.hpp"
#include "../../config/load/cpp_config.hpp"

#include <complex>

using namespace std;

Two_Particle_Interaction::Two_Particle_Interaction(bool complex) {
    if (complex) 
        field_c = (prefix+"_2PI.dat");
    else 
        field_r = (prefix+"_2PI.dat");
    file_found = true;
}

float Two_Particle_Interaction::operator()(Vec k, float w) {
    if (!file_found) return V(k, w);
    if (complex) 
        throw runtime_error("Error: frequency is real, but 2PI is complex");
    return field_r(k, w);
}

complex<float> Two_Particle_Interaction::operator()(Vec k, ::complex<float> w) {
    //if (!file_found) return V(k);
    if (!complex) 
        throw runtime_error("Error: frequency is complex, but 2PI is real");
    return field_c(k, imag(w));
}

