/**
 * Matrix class, implemented in 1D for practicality
 *
 * Author: Griffin Heier
 */

#include "matrix.hpp"
#include "../config/load/cpp_config.hpp"
#include "eigenvec.hpp"
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

Matrix::Matrix() {
    this->size = 0;
}

//Matrix::~Matrix() { delete[] this->vals; }

Matrix::Matrix(int size) {
    this->size = size;
    printv("Matrix size (N): %d\n", size);
    long int size_sqr = static_cast<long>(size) * size;
    printv("Matrix size (NxN): %ld\n", size_sqr);
    float MB = (float)(size_sqr * sizeof(float)) / 1000000;
    if (MB > 100) {
        printv("Allocating %.1lf GB\n", MB / 1024);
    }
    else if (MB > 1) {
        printv("Allocating %.0lf MB\n", MB);
    } else {
        printv("Allocating %.3lf MB\n", MB);
    }
    try {
        vals.resize(size_sqr);
    } catch (const std::bad_alloc& e) {
        std::cerr << "Memory allocation failed: " << e.what() << '\n';
        exit(EXIT_FAILURE);
    }
    printv("Allocated\n");
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            vals[i * size + j] = (i == j) ? 1 : 0;
        }
    }
    printv("Filled\n");
}

Matrix::Matrix(vector<vector<float>> vals) {
    this->size = vals.size();
    vals.resize(size * size);
    for (int i = 0; i < this->size; i++) {
        for (int j = 0; j < this->size; j++) {
            this->vals[i * this->size + j] = vals[i][j];
        }
    }
}

float &Matrix::operator()(int r, int c) {
    return this->vals[r * this->size + c];
}

Matrix Matrix::operator+(const Matrix &k) {
    Matrix result(this->size);
    for (int i = 0; i < this->size; i++) {
        for (int j = 0; j < this->size; j++) {
            result.vals[i * this->size + j] =
                this->vals[i * this->size + j] + k.vals[i * this->size + j];
        }
    }
    return result;
}

Matrix &Matrix::operator+=(const Matrix &k) {
    for (int i = 0; i < this->size; i++) {
        for (int j = 0; j < this->size; j++) {
            this->vals[i * this->size + j] += k.vals[i * this->size + j];
        }
    }
    return *this;
}

Matrix Matrix::operator-(const Matrix &k) {
    Matrix result(this->size);
    for (int i = 0; i < this->size; i++) {
        for (int j = 0; j < this->size; j++) {
            result.vals[i * this->size + j] =
                this->vals[i * this->size + j] - k.vals[i * this->size + j];
        }
    }
    return result;
}

Matrix &Matrix::operator-=(const Matrix &k) {
    for (int i = 0; i < this->size; i++) {
        for (int j = 0; j < this->size; j++) {
            this->vals[i * this->size + j] -= k.vals[i * this->size + j];
        }
    }
    return *this;
}

Matrix Matrix::operator*(const Matrix &k1) {
    if (this->size != k1.size) {
        cout << "Matrix dimensions do not match" << endl;
        exit(1);
        // return Matrix(0);
    }
    Matrix result(this->size);
    for (int i = 0; i < this->size; i++) {
        for (int j = 0; j < this->size; j++) {
            for (int k = 0; k < this->size; k++) {
                result.vals[i * this->size + j] +=
                    this->vals[i * this->size + k] *
                    k1.vals[k * this->size + j];
            }
        }
    }
    return result;
}

Eigenvector Matrix::operator*(Eigenvector &k) {
    Eigenvector result(this->size);
#pragma omp parallel for
    for (int i = 0; i < this->size; i++) {
        for (int j = 0; j < this->size; j++) {
            result.eigenvector[i] +=
                this->vals[i * this->size + j] * k.eigenvector[j];
        }
    }
    return result;
}

Matrix &Matrix::operator*=(float constant) {
    for (int i = 0; i < this->size; i++) {
        for (int j = 0; j < this->size; j++) {
            this->vals[i * this->size + j] *= constant;
        }
    }
    return *this;
}

Matrix Matrix::operator*(float multiple) {
    Matrix result(this->size);
    for (int i = 0; i < this->size; i++) {
        for (int j = 0; j < this->size; j++) {
            result.vals[i * this->size + j] =
                this->vals[i * this->size + j] * multiple;
        }
    }
    return result;
}

Matrix Matrix::operator/(float multiple) {
    Matrix result(this->size);
    for (int i = 0; i < this->size; i++) {
        for (int j = 0; j < this->size; j++) {
            result.vals[i * this->size + j] =
                this->vals[i * this->size + j] / multiple;
        }
    }
    return result;
}

bool Matrix::operator==(const Matrix &k) {
    for (int i = 0; i < this->size; i++) {
        for (int j = 0; j < this->size; j++) {
            if (this->vals[i * this->size + j] != k.vals[i * this->size + j]) {
                return false;
            }
        }
    }
    return true;
}

// std::ostream &operator<<(std::ostream &os, const Matrix &k) {
//     for (int i = 0; i < k.size; i++) {
//         for (int j = 0; j < k.size; j++) {
//             os << k.vals[i * k.size + j] << " ";
//         }
//         os << endl;
//     }
//     return os;
// }
//
// std::istream &operator>>(std::istream &is, Matrix &k) {
//     for (int i = 0; i < k.size; i++) {
//         for (int j = 0; j < k.size; j++) {
//             is >> k.vals[i * k.size + j];
//         }
//     }
//     return is;
// }
