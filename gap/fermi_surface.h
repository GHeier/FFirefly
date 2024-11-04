#pragma once
#ifndef FERMI_SURFACE_H
#define FERMI_SURFACE_H
 
#include <vector>
#include <complex>
#include "vec.h"
#include "surfaces.h"

using namespace std;

struct VecAndEnergy {
    Vec vec;
    float energy;
};

/*
 * @brief Compares two VecAndEnergy structs by energy
 * @param left The first VecAndEnergy struct
 * @param right The second VecAndEnergy struct
 * @return True if the energy of the first struct is less than the energy of the second struct
 */
bool operator<(const VecAndEnergy& left, const VecAndEnergy& right);

/**
 * @brief Tetrahedron method. Defines the surface of a certain energy curve
 *
 * Splits up the surface into a series of cubes, and those into 6 tetrahedrons. The energy surface
 * is approximated to be linear in each tetrahedron, and so it can be calculated this way
 *
 * @param func The function to evaluate
 * @param q The q vector
 * @param s_val The value of the surface
 *
 * @return A vector of Vec structs
 */
vector<Vec> tetrahedron_method(float (*func)(Vec k, Vec q), Vec q, float s_val);
/**
 * @brief Tetrahedron sum. Sums over all the surfaces that are defined continuously, rather than
 * put into a vector
 *
 * This is used for the susceptibility integral, where the surfaces correspond to values of the 
 * denominator, which diverges
 *
 * @param func The function to evaluate
 * @param func_diff The derivative of the function to evaluate
 * @param q The q vector
 * @param svals The values of the surface
 * @param w The frequency
 * @param T The temperature
 *
 * @return The sum of the tetrahedrons
 */
float tetrahedron_sum_continuous(float (*func)(Vec k, Vec q), float (*func_diff)(Vec k, Vec q), Vec q, vector<float> &svals, float w, float T);
/**
 * @brief Analytic Tetrahedron sum, where each tetrahedron is done via an analytic formula
 *
 * This can be accomplished again by approximating the energy surface as linear in each tetrahedron,
 * and then calculating that portion of the integral analytically. Is done at zero temperature, 
 * whereas the method tetrahedron_sum_continuous is not
 *
 * @param q The q vector
 * @param w The frequency
 *
 * @return The sum of the tetrahedrons
 */
float analytic_tetrahedron_sum(Vec q, float w, int num_pts);

#endif
