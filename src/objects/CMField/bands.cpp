#include "bands.hpp"
#include "../../config/load/cpp_config.hpp"
#include "../../hamiltonian/band_structure.hpp"
#include "fields.hpp"
#include "../vec.hpp"
#include "../eigenvec.hpp"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

Bands::Bands() {
    file_found = false;
    nbands = 0;

    string filename = outdir + prefix + "_Hk.h5";
    ifstream file(filename);

    if (file.is_open()) {
        file.close();
        // Load the Hamiltonian as a complex matrix field
        Field_CM H(filename);
        fill_grid(H);
        file_found = true;
        printv("Read in Hamiltonian from %s\n", filename.c_str());
    } else {
        printv("No Hamiltonian file found at %s. Using default band structure\n", filename.c_str());
        nbands = 1;
    }
}

Bands::Bands(Field_CM &H) {
    file_found = true;
    fill_grid(H);
}

void Bands::fill_grid(Field_CM &H) {
    vector<int> nk = H.cmf.data.mesh;
    nbands = H.cmf.data.inds[0];  // First dimension of matrix (number of bands)
    vector<vector<cfloat>> eigs(nbands, vector<cfloat>(nk[0] * nk[1] * nk[2]));
    vector<vector<eigvec>> eigenvectors(nbands, vector<eigvec>(nk[0] * nk[1] * nk[2], eigvec(0)));
    int ind = 0;
    for (int i = 0; i < nk[0]; i++) {
        float x = (float)i / (float)(nk[0] - 1) - 0.5;
        for (int j = 0; j < nk[1]; j++) {
            float y = (float)j / (float)(nk[1] - 1) - 0.5;
            for (int k = 0; k < nk[2]; k++) {
                float z = (float)k / (float)(nk[2] - 1) - 0.5;
                Vec kpoint = brillouin_zone * Vec(x, y, z);
                auto eigvals_and_vecs = H.fulldiag(kpoint, 0.0);
                if (eigvals_and_vecs.empty()) {
                    printf("Error: diagonalization failed at k-point (%f, %f, %f)\n",
                           kpoint.x, kpoint.y, kpoint.z);
                    continue;
                }
                for (int n = 0; n < nbands; n++) {
                    eigs[n][ind] = eigvals_and_vecs[n].eigenvalue;
                    eigenvectors[n][ind] = eigvals_and_vecs[n];
                }
                ind++;
            }
        }
    }
    for (int n = 0; n < nbands; n++) {
        band_fields.push_back(Field_R(eigs[n], nk, H.cmf.data.domain, H.cmf.data.w_points));
        // wavefunctions ... to be implemented when Field_RV is implemented.
    }
}

float Bands::operator()(int n, Vec k) {
    if (n < 1 || n > nbands) {
        printf("Error: band index %d out of range (must be >= 1 and <= %d)\n", n, nbands);
        return 0.0f;
    }
    if (nbands == 0) {
        printf("Error: No bands loaded\n");
        return 0.0f;
    }
    if (!file_found) {
        return epsilon(n, k);
    }

    // Return eigenvalue for band n (n is 1-indexed)
    return band_fields[n - 1](k, 0.0f);
}

float Bands::operator()(Vec k) {
    if (nbands == 0) {
        printf("Error: No bands loaded\n");
        return 0.0f;
    }
    if (!file_found) {
        return epsilon(1, k);
    }

    // Return eigenvalue for band n (n is 1-indexed)
    return band_fields[0](k, 0.0f);
}
