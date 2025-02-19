#include "cmf.hpp"
#include "cmf_types.hpp"

#include "vec.hpp"

complex<float> CMF_CS::operator() (float w) {
    return complex<float>(base(w).real().x, base(w).imag().x);
}

complex<float> CMF_CS::operator() (Vec point, float w) {
    return complex<float>(base(point, w).real().x, base(point, w).imag().x);
}

CMF_CS::CMF_CS() {
    base = CMF();
}

CMF_CS::CMF_CS(CMF base) {
    this->base = base;
}

CMF_CS::CMF_CS(string filename) {
    base = load_CMF_from_file(filename);
}

float CMF_RS::operator() (float w) {
    return base(w).real().x;
}

float CMF_RS::operator() (Vec point, float w) {
    return base(point, w).real().x;
}

CMF_RS::CMF_RS() {
    base = CMF();
}

CMF_RS::CMF_RS(CMF base) {
    this->base = base;
}

CMF_RS::CMF_RS(string filename) {
    base = load_CMF_from_file(filename);
}
