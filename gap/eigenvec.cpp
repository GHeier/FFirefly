<<<<<<< HEAD
/**
 * @file eigenvec.cpp
 *
 * @brief The implementation of the Eigenvector class
 *
 * @author Griffin Heier
 */
=======
>>>>>>> origin/main
#include <string>
#include <vector>
#include <math.h>

<<<<<<< HEAD
#include <bits/stdc++.h>
=======
>>>>>>> origin/main
#include <ctime>
#include <algorithm>

#include "vec.h"
#include "cfg.h"
#include "eigenvec.hpp"
#include "matrix.hpp"

<<<<<<< HEAD
using namespace std;

Eigenvector::Eigenvector() {
    this->size = 0;
    this->eigenvector = new float[0];
    this->eigenvalue = 0;
}

Eigenvector::~Eigenvector() {
    delete[] this->eigenvector;
}

Eigenvector::Eigenvector(int size, bool random) {
    this->size = size;
    this->eigenvector = new float[size];
    if (random) {
        srand(time(0)); 
        for (int i = 0; i < size; i++) {
            this->eigenvector[i] = (float)rand() / RAND_MAX;
        }
    }
}

Eigenvector::Eigenvector(vector<float> eigenvector) {
    this->size = eigenvector.size();
    this->eigenvector = eigenvector.data();
}

Eigenvector::Eigenvector(vector<float> eigenvector, float eigenvalue) {
    this->size = eigenvector.size();
    this->eigenvector = eigenvector.data();
    this->eigenvalue = eigenvalue;
}

float& Eigenvector::operator[](int index) {
=======
Eigenvector::Eigenvector() {
    this->size = 0;
    this->eigenvector = vector<double>(0, 0);
    this->eigenvalue = 0;
}

Eigenvector::Eigenvector(int size) {
    this->size = size;
    this->eigenvector = vector<double>(size);

    std::srand(unsigned(std::time(nullptr)));
    std::generate(this->eigenvector.begin(), this->eigenvector.end(), std::rand);
}

Eigenvector::Eigenvector(vector<double> eigenvector) {
    this->size = eigenvector.size();
    this->eigenvector = eigenvector;
}

Eigenvector::Eigenvector(vector<double> eigenvector, double eigenvalue) {
    this->size = eigenvector.size();
    this->eigenvector = eigenvector;
    this->eigenvalue = eigenvalue;
}

double& Eigenvector::operator[](int index) {
>>>>>>> origin/main
    return this->eigenvector[index];
}

Eigenvector Eigenvector::operator+(const Eigenvector& k) {
    Eigenvector result(this->size);
    for (int i = 0; i < this->size; i++) {
        result.eigenvector[i] = this->eigenvector[i] + k.eigenvector[i];
    }
    return result;
}

Eigenvector Eigenvector::operator-(const Eigenvector& k) {
    Eigenvector result(this->size);
    for (int i = 0; i < this->size; i++) {
        result.eigenvector[i] = this->eigenvector[i] - k.eigenvector[i];
    }
    return result;
}

<<<<<<< HEAD
Eigenvector Eigenvector::operator*(float multiple) {
=======
Eigenvector Eigenvector::operator*(double multiple) {
>>>>>>> origin/main
    Eigenvector result(this->size);
    for (int i = 0; i < this->size; i++) {
        result.eigenvector[i] = this->eigenvector[i] * multiple;
    }
    return result;
}

Eigenvector Eigenvector::operator*(const Eigenvector& k) {
    Eigenvector result(this->size);
    for (int i = 0; i < this->size; i++) {
        result.eigenvector[i] = this->eigenvector[i] * k.eigenvector[i];
    }
    return result;
}

<<<<<<< HEAD
Eigenvector Eigenvector::operator/(float multiple) {
=======
Eigenvector Eigenvector::operator/(double multiple) {
>>>>>>> origin/main
    Eigenvector result(this->size);
    for (int i = 0; i < this->size; i++) {
        result.eigenvector[i] = this->eigenvector[i] / multiple;
    }
    return result;
}

bool Eigenvector::operator==(const Eigenvector& k) {
    for (int i = 0; i < this->size; i++) {
        if (this->eigenvector[i] != k.eigenvector[i]) return false;
    }
    return true;
}

bool operator!=(const Eigenvector& left, const Eigenvector& right) {
    return left.eigenvalue != right.eigenvalue;
}

bool operator<(const Eigenvector& left, const Eigenvector& right) {
    return left.eigenvalue < right.eigenvalue;
}

<<<<<<< HEAD
float Eigenvector::norm() {
    float sum = 0;
=======
double Eigenvector::norm() {
    double sum = 0;
>>>>>>> origin/main
    for (int i = 0; i < this->size; i++) {
        sum += this->eigenvector[i] * this->eigenvector[i];
    }
    return sqrt(sum);
}

void Eigenvector::normalize() {
<<<<<<< HEAD
    float norm = this->norm();
=======
    double norm = this->norm();
>>>>>>> origin/main
    for (int i = 0; i < this->size; i++) {
        this->eigenvector[i] /= norm;
    }
}

<<<<<<< HEAD
float dot(const Eigenvector& left, const Eigenvector& right) {
    float sum = 0;
=======
double dot(const Eigenvector& left, const Eigenvector& right) {
    double sum = 0;
>>>>>>> origin/main
    for (int i = 0; i < left.size; i++) {
        sum += left.eigenvector[i] * right.eigenvector[i];
    }
    return sum;
}

std::ostream& operator<<(std::ostream& os, const Eigenvector& k) {
    os << "[";
    for (int i = 0; i < k.size; i++) {
        os << k.eigenvector[i];
        if (i != k.size - 1) os << ", ";
    }
    os << "]";
    return os;
}

