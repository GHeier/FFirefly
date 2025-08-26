#include "fields.hpp"
#include "cmfield.hpp"

#include "../vec.hpp"

complex<float> Field_C::operator()(int n, float w) {
  return complex<float>(cmf(n, w).real().x, cmf(n, w).imag().x);
}

complex<float> Field_C::operator()(float w) {
  return complex<float>(cmf(w).real().x, cmf(w).imag().x);
}

complex<float> Field_C::operator()(int n, Vec point, float w) {
  return complex<float>(cmf(n, point, w).real().x, cmf(n, point, w).imag().x);
}

complex<float> Field_C::operator()(Vec point, float w) {
  return complex<float>(cmf(point, w).real().x, cmf(point, w).imag().x);
}

Field_C::Field_C() { cmf = CMField(); }

Field_C::Field_C(CMField cmf) { this->cmf = cmf; }

Field_C::Field_C(string filename) { 
    cmf = load_CMField(filename); 
}

float Field_R::operator()(double w) { 
    return cmf(w).real().x; }

float Field_R::operator()(int n, double w) {
    return cmf(n, w).real().x; }

float Field_R::operator()(Vec point, float w) {
    return cmf(point, w).real().x; }

float Field_R::operator()(int n, Vec point, float w) {
  return cmf(n, point, w).real().x;
}

Field_R::Field_R() { 
    cmf = CMField(); 
}

Field_R::Field_R(CMField cmf) { this->cmf = cmf; }

Field_R::Field_R(string filename) { 
    cmf = load_CMField(filename); 
}
