#pragma once
#ifndef POTENTIAL_H_
#define POTENTIAL_H_

#include <vector>
#include "vec.h"
#include "fermi_surface.h"

using namespace std;

struct indices_values {
    vector<int> indices;
    float value;
};

/**
 * @brief Constant potential
 * 
 * @param k1 
 * @param k2 
 * @return float 
 */
float potential_const(Vec k1, Vec k2);

/**
 * @brief Test potential
 * 
 * @param k1 
 * @param k2 
 * @return float 
 */
float potential_test(Vec k1, Vec k2);

/**
 * @brief Phonon Coulomb potential
 * 
 * @param q 
 * @return float 
 */
float phonon_coulomb(Vec q);

/**
 * @brief Scalapino singlet potential (FLEX) done via direct integration
 * 
 * @param k1 
 * @param k2 
 * @param T 
 * @return float 
 */
float potential_scal(Vec k1, Vec k2, float T);

/**
 * @brief Scalapino singlet potential (FLEX) done via chi cube interpolation
 * 
 * @param k1 
 * @param k2 
 * @param w 
 * @param T 
 * @param chi_map 
 * @return float 
 */
float potential_scalapino_cube(Vec k1, Vec k2, float w, float T, const unordered_map<float, vector<vector<vector<float>>>> &chi_map);

/**
 * @brief Scalapino triplet potential (FLEX) done via chi cube interpolation
 * 
 * @param k1 
 * @param k2 
 * @param w 
 * @param T 
 * @param chi_map 
 * @return float 
 */
float potential_scalapino_triplet(Vec k1, Vec k2, float w, float T, const unordered_map<float, vector<vector<vector<float>>>> &chi_map);


#endif
