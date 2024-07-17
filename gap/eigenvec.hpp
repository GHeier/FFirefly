#pragma once
#ifndef eigenvector_class_H
#define eigenvector_class_H

#include <iostream>
#include <string>
#include "vec.h"

<<<<<<< HEAD
/**
 * Eigenvector class
 *
 * This class is used to represent an eigenvector of a matrix.
 * It contains the eigenvector itself, the eigenvalue associated with it, and the size of 
 * the eigenvector.
 */
class Eigenvector {
    public:
        int size;
        float eigenvalue;
        float* eigenvector;

        Eigenvector();
        ~Eigenvector();
        Eigenvector(int size, bool random=false);
        Eigenvector(vector<float> eigenvector);
        Eigenvector(vector<float> eigenvector, float eigenvalue);

        float& operator[](int index);
=======
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
>>>>>>> origin/main

        Eigenvector operator+(const Eigenvector& k);
        Eigenvector operator-(const Eigenvector& k);
        Eigenvector operator*(const Eigenvector& k);
<<<<<<< HEAD
        Eigenvector operator*(float multiple);
        Eigenvector operator/(float multiple);
        bool operator==(const Eigenvector& k);
        bool operator!=(const Eigenvector& k);
        bool operator<(const Eigenvector& k);
        float norm();
=======
        Eigenvector operator*(double multiple);
        Eigenvector operator/(double multiple);
        bool operator==(const Eigenvector& k);
        bool operator!=(const Eigenvector& k);
        bool operator<(const Eigenvector& k);
        double norm();
>>>>>>> origin/main
        void normalize();

};


<<<<<<< HEAD
float dot(const Eigenvector& left, const Eigenvector& right);
=======
double dot(const Eigenvector& left, const Eigenvector& right);
>>>>>>> origin/main
std::ostream& operator<<(std::ostream& os, const Eigenvector& k);

#endif
