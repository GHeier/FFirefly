#pragma once
#ifndef BAND_STRUCTURE_H
#define BAND_STRUCTURE_H

#include <string>
#include "vec.h"

using namespace std;


/**
 * @brief Returns the energy of a free electron in a sphere.
 * 
 * @param k The wave vector.
 * @return float The energy of the electron.
 */
float epsilon_sphere(const Vec k);

/**
 * @brief Returns the Fermi velocity of a free electron in a sphere.
 * 
 * @param k The wave vector.
 * @return Vec The Fermi velocity of the electron.
 */
Vec fermi_velocity_sphere(const Vec k);

/**
 * @brief Returns the energy of an electron in a simple cubic lattice.
 * 
 * @param k The wave vector.
 * @param t The hopping parameter.
 * @param tn The next-nearest neighbor hopping parameter.
 * @return float The energy of the electron.
 */
float epsilon_SC(const Vec k, float t, float tn);

/**
 * @brief Returns the Fermi velocity of an electron in a simple cubic lattice.
 * 
 * @param k The wave vector.
 * @return Vec The Fermi velocity of the electron.
 */
Vec fermi_velocity_SC(const Vec k);

/**
 * @brief Returns the energy of an electron in a layered simple cubic lattice.
 * 
 * @param k The wave vector.
 * @return float The energy of the electron.
 */
float epsilon_SC_layered(const Vec k);

/**
 * @brief Returns the Fermi velocity of an electron in a layered simple cubic lattice.
 * 
 * @param k The wave vector.
 * @return Vec The Fermi velocity of the electron.
 */
Vec fermi_velocity_SC_layered(const Vec k);

#endif
