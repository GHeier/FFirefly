#pragma once
#ifndef SURFACSE_H
#define SURFACSE_H
 
#include <vector>
#include <complex>
#include <functional>

#include "vec.h"

using namespace std;


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
 * @return A vector of Vecs 
 */
vector<Vec> points_from_indices2(function<float(Vec)> func, int i, int j, int k, int divs);

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
vector<Vec> points_in_tetrahedron(function<float(Vec k)> func, float s_val, vector<Vec> points);

/**
 * @brief Checks to see if the energy surface is inside the cube
 *
 * @param s_val The value of the surface
 * @param p The points to check
 *
 * @return True if the surface is inside the cube
 */
bool surface_inside_cube(float s_val, vector<Vec> p);

/**
 * @brief Checks to see if the energy surface is inside the tetrahedron
 *
 * @param s_val The value of the surface
 * @param ep_points The points to check
 *
 * @return True if the surface is inside the tetrahedron
 */
bool surface_inside_tetrahedron(float s_val, vector<Vec> ep_points);

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
vector<Vec> tetrahedron_method2(function<float(Vec k)> func, float s_val);

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

#endif
