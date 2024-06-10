#pragma once
#ifndef eigenvector_class_H
#define eigenvector_class_H

#include <iostream>
#include <string>
#include "vec.h"

class Eigenvector {
    public:
        int size;
        double eigenvalue;
        vector<double> eigenvector;

        Eigenvector();
        Eigenvector(int size);
        Eigenvector(vector<double> eigenvector);
        Eigenvector(vector<double> eigenvector, double eigenvalue);

        double& operator[](int index);

        Eigenvector operator+(const Eigenvector& k);
        Eigenvector operator-(const Eigenvector& k);
        Eigenvector operator*(const Eigenvector& k);
        Eigenvector operator*(double multiple);
        Eigenvector operator/(double multiple);
        bool operator==(const Eigenvector& k);
        bool operator!=(const Eigenvector& k);
        bool operator<(const Eigenvector& k);
        double norm();
        void normalize();

};


double dot(const Eigenvector& left, const Eigenvector& right);
std::ostream& operator<<(std::ostream& os, const Eigenvector& k);

#endif
