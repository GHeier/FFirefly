#include "fields.hpp"

#include <complex>

using namespace std;

class Two_Particle_Interaction {
    public:
        Field_R field_r;
        Field_C field_c;
        bool complex;
        bool file_found;

        Two_Particle_Interaction(bool complex = false);
        ~Two_Particle_Interaction();
        float operator()(Vec k, float w = 0);
        ::complex<float> operator()(Vec k, ::complex<float> w);
};

