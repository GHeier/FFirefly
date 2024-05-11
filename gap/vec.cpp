#include <string>
#include <Eigen/Dense>
#include "vec.h"
#include "cfg.h"

using namespace Eigen;

Vec::Vec(double x, double y, double z, bool is_cartesian, double area, double freq) {
        cartesian = is_cartesian;
        area = area;
        freq = freq;
        //if (dim == 2) z = 0;
        vals(0) = x; vals(1) = y; vals(2) = z;
}

void Vec::to_cartesian() {
    cartesian = true;
    double x = vals(0)*cos(vals(1))*cos(vals(2));
    double y = vals(0)*sin(vals(1))*cos(vals(2));
    double z = vals(0)*sin(vals(2));
    vals(0) = x; vals(1) = y; vals(2) = z;
}
void Vec::to_spherical() {
    cartesian = false;
    double r = vals.norm();
    double theta = atan2(vals(1), vals(0));
    double phi = atan2(vals(2), hypot(vals(0), vals(1)));
    vals(0) = r; vals(1) = theta; vals(2) = phi;
}

Vec string_to_vec(string str) {
    int start = 0;
    int end = str.find(" ");
    Vec result;
    int iter = 0;
    while (end != -1) {
        string temp = str.substr(start, end - start);
        double val = std::stod(temp);
        result.vals(iter) = val;
        start = end + 1;
        end = str.find(" ", start);
        iter++;
    }
    return result;
}

string vec_to_string(Vec k) {
    string hash = "";
    for (int i = 0; i < dim; i++) {
        hash += to_string(std::ceil(k.vals(i)*100.0)/100.0 ) + " ";
    }
    return hash;
}

Vec operator+(const Vec& k, const Vec& q) {
    Vec k_new = k; Vec q_new = q;
    if (k_new.cartesian == false) k_new.to_cartesian();
    if (q_new.cartesian == false) q_new.to_cartesian();
    
    Vec result(k_new.vals(0) + q_new.vals(0), k_new.vals(1) + q_new.vals(1), k_new.vals(2) + q_new.vals(2));

    if (k.cartesian) return result;

    result.to_spherical();
    return result;
}

// Returns in spherical unless both vectors are in cartesian coordinates
Vec operator-(const Vec& k, const Vec& q) {
    Vec k_new = k; Vec q_new = q;
    if (k_new.cartesian == false) k_new.to_cartesian();
    if (q_new.cartesian == false) q_new.to_cartesian();
    
    Vec result(k_new.vals(0) - q_new.vals(0), k_new.vals(1) - q_new.vals(1), k_new.vals(2) - q_new.vals(2));

    if (k.cartesian) return result;

    result.to_spherical();
    return result;
}

Vec operator*(double multiple, const Vec& input) {
    Vec result;
    if (input.cartesian == false) {
        result.vals(0) = input.vals(0) * multiple;
        result.vals(1) = input.vals(1);
        result.vals(2) = input.vals(2);
    }
    else {
        result.vals(0) = input.vals(0) * multiple;
        result.vals(1) = input.vals(1) * multiple;
        result.vals(2) = input.vals(2) * multiple;
    }
    return result;
}

Vec operator*(const Vec& input, double multiple) {
    Vec result;
    if (input.cartesian == false) {
        result.vals(0) = input.vals(0) * multiple;
        result.vals(1) = input.vals(1);
        result.vals(2) = input.vals(2);
    }
    else {
        result.vals(0) = input.vals(0) * multiple;
        result.vals(1) = input.vals(1) * multiple;
        result.vals(2) = input.vals(2) * multiple;
    }
    return result;
}

double operator*(const Vec& left, const Vec& right) {
    double sum = 0;
    for (int i = 0; i < dim; i++) {
        sum += left.vals(i) * right.vals(i);
    }
    return sum;
}

Vec operator/(const Vec& input, double multiple) {
    Vec result;
    if (input.cartesian == false) {
        result.vals(0) = input.vals(0) / multiple;
        result.vals(1) = input.vals(1);
        result.vals(2) = input.vals(2);
    }
    else {
        result.vals(0) = input.vals(0) / multiple;
        result.vals(1) = input.vals(1) / multiple;
        result.vals(2) = input.vals(2) / multiple;
    }
    return result;
}

bool operator==(const Vec& k, const Vec& q) {
    return (fabs( (k-q).vals.norm()) < 0.001);
}

bool operator<(const Vec& p1, const Vec& p2) {
    return p1.area < p2.area;
    return (Vec(2*M_PI,0,0) - p1).vals.norm() < (Vec(2*M_PI,0,0) - p2).vals.norm();
}

std::ostream& operator<<(std::ostream& os, const Vec& k) {
    for (int i = 0; i < dim; i++) {
        os << k.vals(i) << " ";
    }
    return os;
}
