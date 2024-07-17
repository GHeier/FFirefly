/**
 * @file vec.cpp
 *
 * @brief This file contains the implementation of the Vec class
 *
 * @author Griffin Heier
 */

#include <vector>
#include <string>
#include <math.h>
#include "vec.h"
#include "cfg.h"

Vec::Vec(float x, float y, float z, bool is_cartesian, float area, float freq) {
    cartesian = is_cartesian;
    area = area;
    vals[0] = x;
    vals[1] = y;
    vals[2] = z;
    freq = freq;
}

void Vec::to_cartesian() {
    if (!cartesian) {
        float x = vals[0] * sin(vals[1]) * cos(vals[2]);
        float y = vals[0] * sin(vals[1]) * sin(vals[2]);
        float z = vals[0] * cos(vals[1]);
        vals[0] = x;
        vals[1] = y;
        vals[2] = z;
        cartesian = true;
    }
}

void Vec::to_spherical() {
    if (cartesian) {
        float r = sqrt(pow(vals[0], 2) + pow(vals[1], 2) + pow(vals[2], 2));
        float theta = atan2(vals[1], vals[0]);
        float phi = atan2(vals[2], hypot(vals[0], vals[1]));
        vals[0] = r;
        vals[1] = theta;
        vals[2] = phi;
        cartesian = false;
    }
}

Vec string_to_vec(string str) {
    int start = 0;
    int end = str.find(" ");
    Vec result;
    int iter = 0;
    while (end != -1) {
        string temp = str.substr(start, end - start);
        float val = std::stod(temp);
        result.vals[iter] = val;
        start = end + 1;
        end = str.find(" ", start);
        iter++;
    }
    return result;
}

string vec_to_string(Vec k) {
    string result = "";
    for (int i = 0; i < 3; i++) {
        result += to_string(std::ceil(k.vals[i]*100.0)/100.0) + " ";
    }
    return result;
}

Vec operator+(const Vec& k, const Vec& q) {
    Vec k_new = k; Vec q_new = q;
    if (k_new.cartesian == false) k_new.to_cartesian();
    if (q_new.cartesian == false) q_new.to_cartesian();
    
    Vec result(k_new.vals[0] + q_new.vals[0], k_new.vals[1] + q_new.vals[1], k_new.vals[2] + q_new.vals[2]);

    if (k.cartesian) return result;

    result.to_spherical();
    return result;
}

Vec operator-(const Vec& k, const Vec& q) {
    Vec k_new = k; Vec q_new = q;
    if (k_new.cartesian == false) k_new.to_cartesian();
    if (q_new.cartesian == false) q_new.to_cartesian();
    
    Vec result(k_new.vals[0] - q_new.vals[0], k_new.vals[1] - q_new.vals[1], k_new.vals[2] - q_new.vals[2]);

    if (k.cartesian) return result;

    result.to_spherical();
    return result;
}

Vec operator*(const Vec& input, float multiple) {
    Vec result;
    if (input.cartesian == false) {
        result.vals[0] = input.vals[0] * multiple;
        result.vals[1] = input.vals[1];
        result.vals[2] = input.vals[2];
    }
    else {
        result.vals[0] = input.vals[0] * multiple;
        result.vals[1] = input.vals[1] * multiple;
        result.vals[2] = input.vals[2] * multiple;
    }
    return result;
}

Vec operator*(float multiple, const Vec& input) {
    Vec result;
    if (input.cartesian == false) {
        result.vals[0] = input.vals[0] * multiple;
        result.vals[1] = input.vals[1];
        result.vals[2] = input.vals[2];
    }
    else {
        result.vals[0] = input.vals[0] * multiple;
        result.vals[1] = input.vals[1] * multiple;
        result.vals[2] = input.vals[2] * multiple;
    }
    return result;
}

float operator*(const Vec& left, const Vec& right) {
    Vec left_new = left; Vec right_new = right;
    if (left_new.cartesian == false) left_new.to_cartesian();
    if (right_new.cartesian == false) right_new.to_cartesian();
    return left_new.vals[0] * right_new.vals[0] + left_new.vals[1] * right_new.vals[1] + left_new.vals[2] * right_new.vals[2];
}

Vec operator/(const Vec& input, float multiple) {
    Vec result;
    if (input.cartesian == false) {
        result.vals[0] = input.vals[0] / multiple;
        result.vals[1] = input.vals[1];
        result.vals[2] = input.vals[2];
    }
    else {
        result.vals[0] = input.vals[0] / multiple;
        result.vals[1] = input.vals[1] / multiple;
        result.vals[2] = input.vals[2] / multiple;
    }
    return result;
}

float Vec::norm() {
    if (cartesian) {
        return pow(pow(vals[0], 2) + pow(vals[1], 2) + pow(vals[2], 2), 0.5);
    }
    return vals[0];
}

bool operator==(const Vec& k, const Vec& q) {
    Vec k_new = k; Vec q_new = q;
    if (k_new.cartesian == false) k_new.to_cartesian();
    if (q_new.cartesian == false) q_new.to_cartesian();
    return k_new.vals[0] == q_new.vals[0] && k_new.vals[1] == q_new.vals[1] && k_new.vals[2] == q_new.vals[2];
}

bool operator<(const Vec& left, const Vec& right) {
    return left.area <= right.area;
}

std::ostream& operator<<(std::ostream& os, const Vec& k) {
    os << vec_to_string(k);
    return os;
}


