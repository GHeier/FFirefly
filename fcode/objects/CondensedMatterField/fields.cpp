#include "cmf.hpp"
#include "fields.hpp"

#include "../vec.hpp"

complex<float> Field_C::operator() (float w) {
    return complex<float>(base(w).real().x, base(w).imag().x);
}

complex<float> Field_C::operator() (Vec point, float w) {
    printf("point: %f %f %f\n", point.x, point.y, point.z);
    return complex<float>(base(point, w).real().x, base(point, w).imag().x);
}

Field_C::Field_C() {
    base = CMF();
}

Field_C::Field_C(CMF base) {
    this->base = base;
}

Field_C::Field_C(string filename) {
    base = load_CMF_from_file(filename);
}

extern "C" float Field_R::operator() (double w) {
    return base(w).real().x;
}

extern "C" float Field_R::operator() (Vec point, float w) {
    return base(point, w).real().x;
}

Field_R::Field_R() {
    base = CMF();
}

Field_R::Field_R(CMF base) {
    this->base = base;
}

Field_R::Field_R(string filename) {
    base = load_CMF_from_file(filename);
}

extern "C" float test_func(float w) {
    return 2*w;
}
