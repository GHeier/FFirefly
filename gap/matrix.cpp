<<<<<<< HEAD
/**
 * @file matrix.cpp
 *
 * @brief Matrix class, implemented in 1D for practicality
 *
 * @author Griffin Heier
 */

#include <string>
#include <cmath>
=======
#include <string>
>>>>>>> origin/main
#include <vector>
#include "vec.h"
#include "matrix.hpp"
#include "eigenvec.hpp"
#include "cfg.h"

Matrix::Matrix() {
    this->size = 0;
<<<<<<< HEAD
    this->vals = new float[0];
}

Matrix::~Matrix() {
    delete[] this->vals;
=======
    this->vals = vector<vector<double>>(0, vector<double>(0, 0));
>>>>>>> origin/main
}

Matrix::Matrix(int size) {
    this->size = size;
<<<<<<< HEAD
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
=======
    this->vals = vector<vector<double>>(size, vector<double>(size, 0));
    for (int i = 0; i < size; i++) {
        this->vals[i][i] = 1;
    }
}

Matrix::Matrix(vector<vector<double>> vals) {
    this->size = vals.size();
    this->vals = vals;
}

double& Matrix::operator()(int r, int c) {
    return this->vals[r][c];
>>>>>>> origin/main
}

Matrix Matrix::operator+(const Matrix& k) {
    Matrix result(this->size);
    for (int i = 0; i < this->size; i++) {
        for (int j = 0; j < this->size; j++) {
<<<<<<< HEAD
            result.vals[i*this->size + j] = this->vals[i*this->size + j] + k.vals[i*this->size + j];
=======
            result.vals[i][j] = this->vals[i][j] + k.vals[i][j];
>>>>>>> origin/main
        }
    }
    return result;
} 

Matrix Matrix::operator-(const Matrix& k) {
    Matrix result(this->size);
    for (int i = 0; i < this->size; i++) {
        for (int j = 0; j < this->size; j++) {
<<<<<<< HEAD
            result.vals[i*this->size + j] = this->vals[i*this->size + j] - k.vals[i*this->size + j];
=======
            result.vals[i][j] = this->vals[i][j] - k.vals[i][j];
>>>>>>> origin/main
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
<<<<<<< HEAD
                result.vals[i*this->size + j] += this->vals[i*this->size + k] * k1.vals[k*this->size + j];
=======
                result.vals[i][j] += this->vals[i][k] * k1.vals[k][j];
>>>>>>> origin/main
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
<<<<<<< HEAD
            result.eigenvector[i] += this->vals[i*this->size + j] * k.eigenvector[j];
            //if (isnan(result.eigenvector[i])) {
            //    cout << "NAN" << endl;
            //    cout << this->vals[i][j] << " " << k.eigenvector[j] << endl;
            //    exit(1);
            //}
=======
            result.eigenvector[i] += this->vals[i][j] * k.eigenvector[j];
>>>>>>> origin/main
        }
    }
    return result;
}

<<<<<<< HEAD
Matrix Matrix::operator*(float multiple) {
    Matrix result(this->size);
    for (int i = 0; i < this->size; i++) {
        for (int j = 0; j < this->size; j++) {
            result.vals[i*this->size + j] = this->vals[i*this->size + j] * multiple;
=======
Matrix Matrix::operator*(double multiple) {
    Matrix result(this->size);
    for (int i = 0; i < this->size; i++) {
        for (int j = 0; j < this->size; j++) {
            result.vals[i][j] = this->vals[i][j] * multiple;
>>>>>>> origin/main
        }
    }
    return result;
}

<<<<<<< HEAD
Matrix Matrix::operator/(float multiple) {
    Matrix result(this->size);
    for (int i = 0; i < this->size; i++) {
        for (int j = 0; j < this->size; j++) {
            result.vals[i*this->size + j] = this->vals[i*this->size + j] / multiple;
=======
Matrix Matrix::operator/(double multiple) {
    Matrix result(this->size);
    for (int i = 0; i < this->size; i++) {
        for (int j = 0; j < this->size; j++) {
            result.vals[i][j] = this->vals[i][j] / multiple;
>>>>>>> origin/main
        }
    }
    return result;
}

bool Matrix::operator==(const Matrix& k) {
    for (int i = 0; i < this->size; i++) {
        for (int j = 0; j < this->size; j++) {
<<<<<<< HEAD
            if (this->vals[i*this->size + j] != k.vals[i*this->size + j]) {
                return false;
            }
=======
            if (this->vals[i][j] != k.vals[i][j]) return false;
>>>>>>> origin/main
        }
    }
    return true;
}

std::ostream& operator<<(std::ostream& os, const Matrix& k) {
    for (int i = 0; i < k.size; i++) {
        for (int j = 0; j < k.size; j++) {
<<<<<<< HEAD
            os << k.vals[i*k.size + j] << " ";
=======
            os << k.vals[i][j] << " ";
>>>>>>> origin/main
        }
        os << endl;
    }
    return os;
}
