#pragma once
#ifndef FERMI_SURFACE_H
#define FERMI_SURFACE_H
 
#include <vector>
#include "vec.h"

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
 * @brief Calculates the area of a triangle given three points
 *
 * @param k1 The first point
 * @param k2 The second point
 * @param k3 The third point
 *
 * @return The area of the triangle
 */
float triangle_area_from_points(Vec k1, Vec k2, Vec k3);

/**
 * @brief Calculates the points in k-space from a set of indices
 *
 * @param func The function to evaluate
 * @param q The q vector
 * @param i The first index
 * @param j The second index
 * @param k The third index
 *
 * @return A vector of VecAndEnergy structs
 */
vector<VecAndEnergy> points_from_indices(float (*func)(Vec k, Vec q), Vec q, int i, int j, int k, int divs);

/**
 * @brief Picks out the points that are in each tetrahedron
 *
 * @param func The function to evaluate
 * @param q The q vector
 * @param s_val The value of the surface
 * @param points The points to check
 *
 * @return A vector of Vec structs
 */
vector<Vec> points_in_tetrahedron(float (*func)(Vec k, Vec q), Vec q, float s_val, vector<VecAndEnergy> points);

/**
 * @brief Checks to see if the energy surface is inside the cube
 *
 * @param s_val The value of the surface
 * @param p The points to check
 *
 * @return True if the surface is inside the cube
 */
bool surface_inside_cube(float s_val, vector<VecAndEnergy> p);

/**
 * @brief Checks to see if the energy surface is inside the tetrahedron
 *
 * @param s_val The value of the surface
 * @param ep_points The points to check
 *
 * @return True if the surface is inside the tetrahedron
 */
bool surface_inside_tetrahedron(float s_val, vector<VecAndEnergy> ep_points);

/**
 * @brief Calculates the area of the surface in the tetrahedron (which could also be a triangle)
 *
 * @param cp The corner points of the tetrahedron/triangle
 *
 * @return The area of the surface in the tetrahedron/triangle
 */
float area_in_corners(vector<Vec> cp);

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
 * @brief gives index of lowest surface and number of surfaces in cube
 *
 * This is a function designed to find the energy points inside of a surface. When isolating a 
 * single cube, only some surfaces will be contained within that cube. Since performance is
 * essential, we do not want to check every surface, so instead we use a sorted list of energy
 * contours and find the first and last index of the surfaces that are contained within the cube.
 *
 * This is done via a binary search, and the lower index is returned along with the number of 
 * additional surfaces that are contained within the cube.
 *
 * @param L The lower bound of the surface
 * @param U The upper bound of the surface
 * @param sortedList The sorted list of surfaces
 *
 * @return A pair containing the lower index and the number of surfaces contained within the cube
 */
pair<int, int> get_index_and_length(float L, float U, vector<float> &sortedList);

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
float analytic_tetrahedron_sum(Vec q, float w);

#endif
