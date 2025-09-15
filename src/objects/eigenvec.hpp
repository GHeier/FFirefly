#pragma once

#include <iostream>
#include <memory>
#include <vector>

using namespace std;

/**
 * Eigenvector class
 *
 * This class is used to represent an eigenvector of a matrix.
 * It contains the eigenvector itself, the eigenvalue associated with it, and
 * the size of the eigenvector.
 */
class Eigenvector {
  public:
    int size;
    float eigenvalue;
    vector<float> eigenvector;

    Eigenvector();                                   // Default constructor
    Eigenvector(int size, bool random = false);              // Parameterized constructor

    float& operator[](int index);          

    Eigenvector &operator+=(const Eigenvector &k);
    Eigenvector &operator-=(const Eigenvector &k);
    Eigenvector operator+(const Eigenvector &k);
    Eigenvector operator-(const Eigenvector &k);
    Eigenvector &operator*=(float multiple);
    Eigenvector &operator/=(float multiple);
    Eigenvector operator*(const Eigenvector &k);
    Eigenvector operator*(float multiple);
    Eigenvector operator/(float multiple);
    bool operator==(const Eigenvector &k);
    bool operator!=(const Eigenvector &k);
    bool operator<(const Eigenvector &k);
    float norm();
    void normalize();
};

bool descending_eigenvalues(const Eigenvector &left, const Eigenvector &right);
float dot(const Eigenvector &left, const Eigenvector &right);
// std::ostream& operator<<(std::ostream& os, const Eigenvector& k);
