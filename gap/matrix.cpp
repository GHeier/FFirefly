#include <string>
#include <vector>
#include "vec.h"
#include "matrix.hpp"
#include "eigenvec.hpp"
#include "cfg.h"

Matrix::Matrix() {
    this->size = 0;
    this->vals = vector<vector<double>>(0, vector<double>(0, 0));
}

Matrix::Matrix(int size) {
    this->size = size;
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
}

Matrix Matrix::operator+(const Matrix& k) {
    Matrix result(this->size);
    for (int i = 0; i < this->size; i++) {
        for (int j = 0; j < this->size; j++) {
            result.vals[i][j] = this->vals[i][j] + k.vals[i][j];
        }
    }
    return result;
} 

Matrix Matrix::operator-(const Matrix& k) {
    Matrix result(this->size);
    for (int i = 0; i < this->size; i++) {
        for (int j = 0; j < this->size; j++) {
            result.vals[i][j] = this->vals[i][j] - k.vals[i][j];
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
                result.vals[i][j] += this->vals[i][k] * k1.vals[k][j];
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
            result.eigenvector[i] += this->vals[i][j] * k.eigenvector[j];
        }
    }
    return result;
}

Matrix Matrix::operator*(double multiple) {
    Matrix result(this->size);
    for (int i = 0; i < this->size; i++) {
        for (int j = 0; j < this->size; j++) {
            result.vals[i][j] = this->vals[i][j] * multiple;
        }
    }
    return result;
}

Matrix Matrix::operator/(double multiple) {
    Matrix result(this->size);
    for (int i = 0; i < this->size; i++) {
        for (int j = 0; j < this->size; j++) {
            result.vals[i][j] = this->vals[i][j] / multiple;
        }
    }
    return result;
}

bool Matrix::operator==(const Matrix& k) {
    for (int i = 0; i < this->size; i++) {
        for (int j = 0; j < this->size; j++) {
            if (this->vals[i][j] != k.vals[i][j]) return false;
        }
    }
    return true;
}

std::ostream& operator<<(std::ostream& os, const Matrix& k) {
    for (int i = 0; i < k.size; i++) {
        for (int j = 0; j < k.size; j++) {
            os << k.vals[i][j] << " ";
        }
        os << endl;
    }
    return os;
}
