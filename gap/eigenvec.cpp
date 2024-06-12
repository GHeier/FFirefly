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
    this->eigenvector = vector<double>(0, 0);
    this->eigenvalue = 0;
}

Eigenvector::Eigenvector(int size, bool random) {
    this->size = size;
    this->eigenvector = vector<double>(size, 0);
    if (random) {
        srand(time(0)); 
        for (int i = 0; i < size; i++) {
            this->eigenvector[i] = (double)rand() / RAND_MAX;
        }
    }
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

Eigenvector Eigenvector::operator*(double multiple) {
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

Eigenvector Eigenvector::operator/(double multiple) {
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

double Eigenvector::norm() {
    double sum = 0;
    for (int i = 0; i < this->size; i++) {
        sum += this->eigenvector[i] * this->eigenvector[i];
    }
    return sqrt(sum);
}

void Eigenvector::normalize() {
    double norm = this->norm();
    for (int i = 0; i < this->size; i++) {
        this->eigenvector[i] /= norm;
    }
}

double dot(const Eigenvector& left, const Eigenvector& right) {
    double sum = 0;
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

