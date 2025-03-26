#pragma once

#include "../cmdata.hpp"

bool CMData_tests();
complex<Vec> func_d(Vec x, int dim = 3);
float indf(int i);

bool readwrite_1d();
bool readwrite_2d();
bool readwrite_3d();
bool readwrite_1d_complex();
bool readwrite_2d_complex();
bool readwrite_3d_complex();
bool readwrite_0d_with_w();
bool readwrite_1d_with_w();
bool readwrite_2d_with_w();
bool readwrite_3d_with_w();
bool readwrite_1d_vector();
bool readwrite_2d_vector();
bool readwrite_3d_vector();
bool readwrite_1d_complexvector();
bool readwrite_2d_complexvector();
bool readwrite_3d_complexvector();
bool readwrite_0d_complexvector_with_w();
bool readwrite_1d_complexvector_with_w();
bool readwrite_2d_complexvector_with_w();
bool readwrite_3d_complexvector_with_w();
bool readwrite_1d_with_n();
bool readwrite_2d_with_n();
bool readwrite_3d_with_n();
