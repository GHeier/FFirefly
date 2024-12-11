/**
 * The implementation of the Eigenvector class
 *
 * Author: Griffin Heier
 */

#include <math.h>

#include <bits/stdc++.h>
#include <ctime>
#include <memory>

#include "../config/load/cpp_config.hpp"
#include "eigenvec.hpp"

using namespace std;

Eigenvector::Eigenvector() {
    this->size = 0;
    this->eigenvector = nullptr;
    this->eigenvalue = 0;
}

Eigenvector::~Eigenvector() {}

Eigenvector::Eigenvector(int size, bool random) {
    this->size = size;
    eigenvector = make_unique<float[]>(size);
    if (random) {
        srand(time(0)); 
        for (int i = 0; i < size; i++) {
            this->eigenvector[i] = (float)rand() / RAND_MAX;
        }
    }
}

Eigenvector::Eigenvector(const Eigenvector& other) {
    this->size = other.size;
    this->eigenvalue = other.eigenvalue;
    eigenvector = make_unique<float[]>(size);
    for (int i = 0; i < size; i++) {
        this->eigenvector[i] = other.eigenvector[i];
    }
}

Eigenvector Eigenvector::operator=(const Eigenvector& other) {
    if (this != &other) {
        size = other.size;
        eigenvalue = other.eigenvalue;
        eigenvector = make_unique<float[]>(size);
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

bool Eigenvector::operator!=(const Eigenvector& k) {
    return !(*this == k);
}

bool Eigenvector::operator<(const Eigenvector& k) {
    return this->eigenvalue < k.eigenvalue;
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

bool descending_eigenvalues(const Eigenvector& left, const Eigenvector& right) {
    return left.eigenvalue > right.eigenvalue;
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

