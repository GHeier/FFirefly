#include "cmf.hpp"
#include "fields.hpp"

#include "../vec.hpp"

complex<float> Field_C::operator() (float w) {
    return complex<float>(cmf(w).real().x, cmf(w).imag().x);
}

complex<float> Field_C::operator() (Vec point, float w) {
    printf("point: %f %f %f\n", point.x, point.y, point.z);
    return complex<float>(cmf(point, w).real().x, cmf(point, w).imag().x);
}

Field_C::Field_C() {
    cmf = CMF();
}

Field_C::Field_C(CMF cmf) {
    this->cmf = cmf;
}

Field_C::Field_C(string filename) {
    cmf = load_CMF_from_file(filename);
}

extern "C" float Field_R::operator() (double w) {
    return cmf(w).real().x;
}

extern "C" float Field_R::operator() (Vec point, float w) {
    return cmf(point, w).real().x;
}

Field_R::Field_R() {
    cmf = CMF();
}

Field_R::Field_R(CMF cmf) {
    this->cmf = cmf;
}

Field_R::Field_R(string filename) {
    cmf = load_CMF_from_file(filename);
}

extern "C" float test_func(float w) {
    return 2*w;
}
