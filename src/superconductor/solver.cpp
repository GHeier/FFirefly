#include <ctime>
#include <iostream>
#include <unistd.h>

#include <boost/functional/hash.hpp>
#include <math.h>
#include <openblas/lapacke.h>
#include <unordered_map>

// #include <lambda_lanczos/lambda_lanczos.hpp>
#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/tools/roots.hpp>

#include "../algorithms/linear_algebra.hpp"
#include "../config/load/cpp_config.hpp"
#include "../hamiltonian/band_structure.hpp"
#include "../hamiltonian/interaction.hpp"
#include "../objects/eigenvec.hpp"
#include "../objects/matrix.hpp"
#include "../objects/vec.hpp"
#include "cfg.hpp"
#include "frequency_inclusion.hpp"
#include "matrix_creation.hpp"

/*
 * f_singlet is the part of the linearized BCS gap equation:
 *  epsilon * Delta_k = sum { V_kk' Delta_k' tanh(E_k / 2*T) / (2*E_k) }
 *  epsilon * Delta_k = sum { V_kk' Delta_k' f_singlet(E_k, T) }
 *  f_singlet(E_k, T) = tanh(E_k / 2*T) / (2*E_k)
 */
float f_singlet(float x, float T) {
    if (abs(x) < 0.00001)
        return 1 / (4 * T);
    return tanh(x / (2 * T)) / (2 * x);
}

// Integral of f_singlet over E_k from -wD to wD
// wD is the debye frequency
float f_singlet_integral(float T) {
    auto f = [T](float x) { return f_singlet(x, T); };
    float integral =
        2 * boost::math::quadrature::gauss<float, 7>::integrate(f, 0, wc);
    // printf("Singlet Integral: %.5f\n", integral);
    return integral;
}

float get_Tc_FS_only(double eig) {
    float lower = 0.000001f; // Define appropriate lower bound
    float upper = 1.0f;      // Define appropriate upper bound
    // printf("T=0.00001 : f = %f\n", f_singlet_integral(0.00001));
    // printf("T=1.0 : f = %f\n", f_singlet_integral(1.0));

    // double T_ex = 0.025;
    // double f_ex = f_singlet_integral(T_ex);
    // printf("Expected: %f, Calculated %f\n", 1.512, f_ex);

    auto result = boost::math::tools::bisect(
        [eig](float T) {
            return f_singlet_integral(T) * eig - 1.0f;
        }, // root function
        lower, upper,
        [](float l, float u) {
            return (u - l) < 0.000001f;
        } // termination condition
    );

    float root = (result.first + result.second) / 2.0f;
    return root;
}

// Returns the highest eigenvalue-1 of a given matrix V at temperature T
// This is the function used for root finding by get_Tc
// By finding the root of eig-1, we find the temperature where eig=1
float f(vector<Vec> k, float T,
        const unordered_map<float, vector<vector<vector<float>>>> &cube_map) {
    cout << "\nTemperature point: " << T << endl;
    Matrix P(k.size());
    create_P(P, k);
    float f_integrated = f_singlet_integral(T);

    vector<Eigenvector> answers = power_iteration(P, 0.0001);
    float eig = answers[answers.size() - 1].eigenvalue;
    eig *= f_integrated;

    cout << "Calculated Eigenvalue: " << eig << endl;
    return eig - 1;
}

// Returns the temperature where eig=1
// Uses the f() function above to achieve that, just finds the root of the
// function
float get_Tc(
    vector<Vec> k,
    const unordered_map<float, vector<vector<vector<float>>>> &cube_map) {
    float lower = 0.0005;
    float upper = 1;

    cout << "Determining if Tc exists...\n";
    float max_eig = f(k, lower, cube_map);
    cout << "Maximum eigenvalue is: " << max_eig + 1 << endl;
    if (max_eig > 0) {
        cout << "Tc is less than 5K. Returning 0\n";
        exit(1);
    }
    cout << "Tc exists. Calculating exact Critical Temperature...\n";

    auto x = boost::math::tools::bisect(
        [k, cube_map](float T) { return f(k, T, cube_map); }, lower, upper,
        [=](float lower, float upper) { return upper - lower < 0.0001; });

    cout << "Lower: " << x.first << " Upper: " << x.second << endl;
    float root = (x.second + x.first) / 2;
    return root;
}

float coupling_calc(vector<Vec> &FS, float T) {
    cout << "Calculating Coupling Constant...\n";
    int size = FS.size();
    float DOS = get_DOS(FS);
    // auto cube_map = chi_cube_freq(T, mu);
    // auto cube = chi_cube(T, mu, DOS, 0);
    float f_integrated = f_singlet_integral(T);

    float lambda = 0, normalization = 0;
    auto wave = [](Vec k) { return cos(k(1)) - cos(k(0)); };
    for (int i = 0; i < size; i++) {
        Vec k1 = FS[i];
        for (int j = 0; j < size; j++) {
            Vec k2 = FS[j];
            // lambda += - k1.area * k2.area / vp(k1.n, k1) * wave(k1) * V(k1,
            // k2, T, cube) / vp(k2.n, k2) * wave(k2);
            lambda += V(k1 - k2) * k1.area * k2.area;
        }
        normalization += pow(wave(k1), 2) * k1.area / vp(k1.n, k1);
    }
    cout << "Normalization: " << normalization << endl;
    cout << "Lambda: " << lambda << endl;
    return lambda / normalization * (2 / pow(2 * M_PI, 3));
}
