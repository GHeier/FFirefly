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
#include "../cfg.h"

Vec::Vec(float _x, float _y, float _z, float _w, float _area, int _dimension) {
    dimension = _dimension;
    area = _area;
    x = _x;
    y = _y;
    z = _z;
    w = _w;
}

float& Vec::operator()(int i) {
    if (i == 0) return x;
    if (i == 1) return y;
    if (i == 2) return z;
    if (i == 3) return z;
    printf("Invalid index for Vec\n");
    exit(1);
}


Vec string_to_vec(string str) {
    int start = 0;
    int end = str.find(" ");
    Vec result;
    int iter = 0;
    while (end != -1) {
        string temp = str.substr(start, end - start);
        float val = std::stod(temp);
        if (iter == 0) result.x = val;
        if (iter == 1) result.y = val;
        if (iter == 2) result.z = val;
        if (iter == 3) result.w = val;
        start = end + 1;
        end = str.find(" ", start);
        iter++;
    }
    return result;
}

vector<float> unpack_string(string str) {
    vector<string> result(4);
    int start = 0;
    int end = str.find(" ");
    int iter = 0;
    while (end != -1) {
        string temp = str.substr(start, end - start);
        result[iter] = temp;
        start = end + 1;
        end = str.find(" ", start);
        iter++;
    }
    result[3] = str.substr(start, str.length() - start);
    vector<float> result_float(4);
    for (int i = 0; i < 4; i++) {
        if (result[i] == "") {
            result_float[i] = 0;
            continue;
        }
        result_float[i] = std::stod(result[i]);
    }
    return result_float;
}

string vec_to_string(Vec k) {
    string result = to_string(k.x) + " " 
        + to_string(k.y) + " " 
        + to_string(k.z) + " " 
        + to_string(k.w);
    return result;
}

Vec operator+(const Vec& k, const Vec& q) {
    Vec result(k.x + q.x, k.y + q.y, k.z + q.z, k.w + q.w);
    return result;
}

Vec operator-(const Vec& k, const Vec& q) {
    Vec result(k.x - q.x, k.y - q.y, k.z - q.z, k.w - q.w);
    return result;
}

Vec operator*(const Vec& input, float multiple) {
    Vec result(input.x * multiple, input.y * multiple, input.z * multiple, input.w * multiple);
    return result;
}

Vec operator*(float multiple, const Vec& input) {
    Vec result(input.x * multiple, input.y * multiple, input.z * multiple, input.w * multiple);
    return result;
}

float operator*(const Vec& left, const Vec& right) {
    Vec left_new = left; Vec right_new = right;
    return left_new.x * right_new.x + left_new.y * right_new.y + left_new.z * right_new.z + left_new.w * right_new.w;
}

Vec operator/(const Vec& input, float multiple) {
    Vec result(input.x / multiple, input.y / multiple, input.z / multiple, input.w / multiple);
    return result;
}

float Vec::norm() {
    if (dimension == 1) return sqrt(x * x);
    if (dimension == 2) return sqrt(x * x + y * y);
    if (dimension == 3) return sqrt(x * x + y * y + z * z);
    return sqrt(x * x + y * y + z * z + w * w);
}

bool operator==(const Vec& k, const Vec& q) {
    return k.x == q.x && k.y == q.y && k.z == q.z && k.w == q.w;
}

bool operator<(const Vec& left, const Vec& right) {
    return left.w < right.w;
}

std::ostream& operator<<(std::ostream& os, const Vec& k) {
    if (k.dimension == 1) os << k.x;
    if (k.dimension == 2) os << k.x << " " << k.y;
    if (k.dimension == 3) os << k.x << " " << k.y << " " << k.z;
    if (k.dimension == 4) os << k.x << " " << k.y << " " << k.z << " " << k.w;
    return os;
}


