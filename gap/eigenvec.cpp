/**
 * @file eigenvec.cpp
 *
 * @brief The implementation of the Eigenvector class
 *
 * @author Griffin Heier
 */

#include <string>
#include <vector>
#include <math.h>

#include <bits/stdc++.h>
#include <ctime>
#include <algorithm>

#include "vec.h"
#include "cfg.h"
#include "eigenvec.hpp"
#include "matrix.hpp"

using namespace std;

Eigenvector::Eigenvector() {
    this->size = 0;
    this->eigenvector = nullptr;
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

Eigenvector::Eigenvector(const Eigenvector& other) {
    this->size = other.size;
    this->eigenvalue = other.eigenvalue;
    this->eigenvector = new float[size];
    for (int i = 0; i < size; i++) {
        this->eigenvector[i] = other.eigenvector[i];
    }
}

Eigenvector Eigenvector::operator=(const Eigenvector& other) {
    if (this != &other) {
        delete[] eigenvector;
        size = other.size;
        eigenvalue = other.eigenvalue;
        eigenvector = new float[size];
        for (int i = 0; i < size; i++) {
            this->eigenvector[i] = other.eigenvector[i];
        }
    }
    return *this;
}

float& Eigenvector::operator[](int index) {
    return eigenvector[index];
}

Eigenvector& Eigenvector::operator+=(const Eigenvector& k) {
    for (int i = 0; i < this->size; i++) {
        this->eigenvector[i] += k.eigenvector[i];
    }
    return *this;
}

Eigenvector& Eigenvector::operator-=(const Eigenvector& k) {
    for (int i = 0; i < this->size; i++) {
        this->eigenvector[i] -= k.eigenvector[i];
    }
    return *this;
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

Eigenvector& Eigenvector::operator*=(float constant) {
    for (int i = 0; i < this->size; i++) {
        this->eigenvector[i] *= constant;
    }
    return *this;
}

Eigenvector& Eigenvector::operator/=(float constant) {
    for (int i = 0; i < this->size; i++) {
        this->eigenvector[i] /= constant;
    }
    return *this;
}

Eigenvector Eigenvector::operator*(float multiple) {
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

Eigenvector Eigenvector::operator/(float multiple) {
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

float Eigenvector::norm() {
    float sum = 0;
    for (int i = 0; i < this->size; i++) {
        sum += this->eigenvector[i] * this->eigenvector[i];
    }
    return sqrt(sum);
}

void Eigenvector::normalize() {
    float norm = this->norm();
    for (int i = 0; i < this->size; i++) {
        this->eigenvector[i] /= norm;
    }
}

float dot(const Eigenvector& left, const Eigenvector& right) {
    float sum = 0;
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

