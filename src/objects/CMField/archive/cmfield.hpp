/*
 * This file has been deprecated. CMField class has been archived.
 * The interpolation functions (CMF_search_*) have been moved to field_funcs.hpp
 * Use Field_C, Field_R, Field_CM, or Field_RM instead of CMField.
 *
 * For legacy code, see archive/old_cmfield.hpp
 *
 * Author: Griffin Heier
 */
#pragma once

#include <string>
#include <vector>

using namespace std;

// For backward compatibility, include the field function declarations
#include "field_funcs.hpp"

// Legacy save function still used by other parts of codebase
void save_to_field(string filename, vector<vector<vector<float>>> &values, vector<vector<float>> &domain, vector<int> mesh, vector<float> w_points, bool is_complex, bool is_vector);
