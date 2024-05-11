#pragma once
#ifndef vec_class_H
#define vec_class_H

#include <iostream>
#include <string>
#include <Eigen/Dense>

using namespace std;

class Vec {
    public:
        Eigen::Vector3d vals = Eigen::Vector3d(0,0,0);
        bool cartesian = false;
        double area = 0;
        double freq = 0;

        Vec(double x = 0, double y = 0, double z = 0, bool is_cartesian=true, double area = 0, double freq = 0);
        void to_cartesian();
        void to_spherical();

};

Vec string_to_vec(string str);
string vec_to_string(Vec k);
Vec operator+(const Vec& k, const Vec& q);
Vec operator-(const Vec& k, const Vec& q);
Vec operator*(double multiple, const Vec &input);
Vec operator*(const Vec &input, double multiple);
double operator*(const Vec& left, const Vec& right);
Vec operator/(const Vec &input, double multiple);
bool operator==(const Vec& k, const Vec& q);
bool operator<(const Vec& left, const Vec& right);
std::ostream& operator<<(std::ostream& os, const Vec& k);

#endif
