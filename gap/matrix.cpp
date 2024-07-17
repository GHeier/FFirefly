/**
 * @file matrix.cpp
 *
 * @brief Matrix class, implemented in 1D for practicality
 *
 * @author Griffin Heier
 */

#include <string>
#include <cmath>
#include <vector>
#include "vec.h"
#include "matrix.hpp"
#include "eigenvec.hpp"
#include "cfg.h"

Matrix::Matrix() {
    this->size = 0;
    this->vals = new float[0];
}

Matrix::~Matrix() {
    delete[] this->vals;
}

Matrix::Matrix(int size) {
    this->size = size;
    this->vals = new float[size*size]; 
    for (int i = 0; i < size; i++) {
        this->vals[i*size + i] = 1;
    }
}

Matrix::Matrix(vector<vector<float>> vals) {
    this->size = vals.size();
    this->vals = new float[this->size*this->size];
    for (int i = 0; i < this->size; i++) {
        for (int j = 0; j < this->size; j++) {
            this->vals[i*this->size + j] = vals[i][j];
        }
    }
}

float& Matrix::operator()(int r, int c) {
    return this->vals[r*this->size + c];
}

Matrix Matrix::operator+(const Matrix& k) {
    Matrix result(this->size);
    for (int i = 0; i < this->size; i++) {
        for (int j = 0; j < this->size; j++) {
            result.vals[i*this->size + j] = this->vals[i*this->size + j] + k.vals[i*this->size + j];
        }
    }
    return result;
} 

Matrix Matrix::operator-(const Matrix& k) {
    Matrix result(this->size);
    for (int i = 0; i < this->size; i++) {
        for (int j = 0; j < this->size; j++) {
            result.vals[i*this->size + j] = this->vals[i*this->size + j] - k.vals[i*this->size + j];
        }
    }
    return result;
}

Matrix Matrix::operator*(const Matrix& k1) {
    if (this->size != k1.size) {
        cout << "Matrix dimensions do not match" << endl;
        exit(1);
        //return Matrix(0);
    }
    Matrix result(this->size);
    for (int i = 0; i < this->size; i++) {
        for (int j = 0; j < this->size; j++) {
            for (int k = 0; k < this->size; k++) {
                result.vals[i*this->size + j] += this->vals[i*this->size + k] * k1.vals[k*this->size + j];
            }
        }
    }
    return result;
}

Eigenvector Matrix::operator*(Eigenvector& k) {
    Eigenvector result(this->size);
    #pragma omp parallel for
    for (int i = 0; i < this->size; i++) {
        for (int j = 0; j < this->size; j++) {
            result.eigenvector[i] += this->vals[i*this->size + j] * k.eigenvector[j];
            //if (isnan(result.eigenvector[i])) {
            //    cout << "NAN" << endl;
            //    cout << this->vals[i][j] << " " << k.eigenvector[j] << endl;
            //    exit(1);
            //}
        }
    }
    return result;
}

Matrix Matrix::operator*(float multiple) {
    Matrix result(this->size);
    for (int i = 0; i < this->size; i++) {
        for (int j = 0; j < this->size; j++) {
            result.vals[i*this->size + j] = this->vals[i*this->size + j] * multiple;
        }
    }
    return result;
}

Matrix Matrix::operator/(float multiple) {
    Matrix result(this->size);
    for (int i = 0; i < this->size; i++) {
        for (int j = 0; j < this->size; j++) {
            result.vals[i*this->size + j] = this->vals[i*this->size + j] / multiple;
        }
    }
    return result;
}

bool Matrix::operator==(const Matrix& k) {
    for (int i = 0; i < this->size; i++) {
        for (int j = 0; j < this->size; j++) {
            if (this->vals[i*this->size + j] != k.vals[i*this->size + j]) {
                return false;
            }
        }
    }
    return true;
}

std::ostream& operator<<(std::ostream& os, const Matrix& k) {
    for (int i = 0; i < k.size; i++) {
        for (int j = 0; j < k.size; j++) {
            os << k.vals[i*k.size + j] << " ";
        }
        os << endl;
    }
    return os;
}
