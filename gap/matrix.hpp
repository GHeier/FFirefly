#pragma once
#ifndef matrix_class_H
#define matrix_class_H

#include <iostream>
#include <string>
#include "vec.h"
#include "eigenvec.hpp"

using namespace std;

class Matrix {
    public:
        int size;
        vector<vector<double>> vals;
        Matrix();
        Matrix(int size);
        Matrix(vector<vector<double>> vals);

        double& operator()(int r, int c);

        Matrix operator+(const Matrix& k);
        Matrix operator-(const Matrix& k);
        Matrix operator*(const Matrix& k);
        Eigenvector operator*(Eigenvector& k);
        Matrix operator*(double multiple);
        Matrix operator/(double multiple);
        bool operator==(const Matrix& k);
};

std::ostream& operator<<(std::ostream& os, const Matrix& k);

#endif
