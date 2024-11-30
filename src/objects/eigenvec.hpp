#pragma once
#ifndef eigenvector_class_H
#define eigenvector_class_H

#include <iostream>
#include <memory>
#include <vector>

using namespace std;

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
        unique_ptr<float[]> eigenvector;

        Eigenvector();
        ~Eigenvector();
        Eigenvector(int size, bool random=false);
        Eigenvector(vector<float> eigenvector);
        Eigenvector(vector<float> eigenvector, float eigenvalue);

        Eigenvector(const Eigenvector& other);
        Eigenvector operator=(const Eigenvector& other);
        float& operator[](int index);

        Eigenvector& operator+=(const Eigenvector& k);
        Eigenvector& operator-=(const Eigenvector& k);
        Eigenvector operator+(const Eigenvector& k);
        Eigenvector operator-(const Eigenvector& k);
        Eigenvector& operator*=(float multiple);
        Eigenvector& operator/=(float multiple);
        Eigenvector operator*(const Eigenvector& k);
        Eigenvector operator*(float multiple);
        Eigenvector operator/(float multiple);
        bool operator==(const Eigenvector& k);
        bool operator!=(const Eigenvector& k);
        bool operator<(const Eigenvector& k);
        float norm();
        void normalize();

};


bool descending_eigenvalues(const Eigenvector& left, const Eigenvector& right);
float dot(const Eigenvector& left, const Eigenvector& right);
std::ostream& operator<<(std::ostream& os, const Eigenvector& k);

#endif
