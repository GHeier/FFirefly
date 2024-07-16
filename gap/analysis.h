#pragma once
#ifndef ANALYSIS_H_
#define ANALYSIS_H_

#include "cfg.h"
#include "fermi_surface.h"
#include "potential.h"
#include "vec.h"
#include "save_data.h"
#include "utilities.h"

using namespace std;

/**
 * @brief Plot the hopping integral ratios.
 *
 * The larger the nearest neighbor hopping, generally the farther shifted the Fermi Surface
 * nesting point is, which leads to a stronger Tc. This project was not continued
 *
 */
void plot_hopping_integral_ratios();


#endif
