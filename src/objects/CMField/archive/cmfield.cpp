/*
 * This file has been deprecated. CMField class has been archived.
 * The interpolation functions (CMF_search_*) have been moved to field_funcs.cpp
 * Use Field_C, Field_R, Field_CM, or Field_RM instead of CMField.
 *
 * For legacy code, see archive/old_cmfield.cpp
 *
 * Author: Griffin Heier
 */

#include "field_funcs.hpp"
#include "base_data.hpp"
#include "../vec.hpp"
#include <string>
#include <vector>
#include <complex>

using namespace std;
using cfloat = complex<float>;

// Legacy save function still used by other parts of codebase
void save_to_field(string filename, vector<vector<vector<float>>> &values, vector<vector<float>> &domain, vector<int> mesh, vector<float> w_points, bool is_complex, bool is_vector) {
    int nbnd = values.size();
    int nkpt = values[0].size();

    // Convert vector<vector<vector<float>>> to BaseData::DataVariant format
    // values[band][kpoint][data] where data contains real (and imaginary if complex)

    vector<cfloat> flat_data;

    if (nbnd == 1) {
        // Single band case - use vector<cfloat>
        for (int k = 0; k < nkpt; k++) {
            if (is_complex) {
                // values[0][k] = [real, imag]
                cfloat val(values[0][k][0], values[0][k][1]);
                flat_data.push_back(val);
            } else {
                // values[0][k] = [real]
                cfloat val(values[0][k][0], 0.0f);
                flat_data.push_back(val);
            }
        }
        BaseData::DataVariant data_variant = flat_data;
        save_data(filename, data_variant, is_complex, mesh, domain, w_points, 0, 1);
    } else {
        // Multi-band case - treat bands as indices
        vector<vector<cfloat>> indexed_data;

        for (int k = 0; k < nkpt; k++) {
            vector<cfloat> band_vals(nbnd);
            for (int b = 0; b < nbnd; b++) {
                if (is_complex) {
                    band_vals[b] = cfloat(values[b][k][0], values[b][k][1]);
                } else {
                    band_vals[b] = cfloat(values[b][k][0], 0.0f);
                }
            }
            indexed_data.push_back(band_vals);
        }

        BaseData::DataVariant data_variant = indexed_data;
        save_data(filename, data_variant, is_complex, mesh, domain, w_points, 1, nbnd);
    }
}
