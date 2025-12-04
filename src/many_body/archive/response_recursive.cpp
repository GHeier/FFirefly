#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <array>

using namespace std;

// Global parameters
double t = 1.0;
double mu = -1.0;
const double EPSILON = 1e-10;

struct Vec {
    double x, y;
    Vec() : x(0), y(0) {}
    Vec(double x_, double y_) : x(x_), y(y_) {}
    Vec operator+(const Vec& v) const { return Vec(x + v.x, y + v.y); }
    Vec operator-(const Vec& v) const { return Vec(x - v.x, y + v.y); }
};

Vec get_k(int i, int j, int Nk) {
    return Vec((2.0 * M_PI / Nk) * i, (2.0 * M_PI / Nk) * j);
}

double epsilon(Vec k) {
    return -2.0 * t * (cos(k.x) + cos(k.y));
}

double fermi_dirac(double e) {
    return (e < mu) ? 1.0 : 0.0;
}

// Triangle integral formula from Appendix C of the recursive paper
// For W(k) = 1/D(k) with linearized D inside triangle
double phi_function(double a, double b, double c) {
    // phi(a,b,c) = -b/(2(a-b)(c-b)) - b^2*(-log|a|+log|b|)/(2(a-b)^2(c-b))

    // Handle limiting cases
    if (abs(a - b) < EPSILON) {  // a -> b
        if (abs(b - c) < EPSILON) return 0.0;
        return 1.0 / (4.0 * b - 4.0 * c);
    }

    if (abs(b - c) < EPSILON) return 0.0;  // b -> c handled by caller

    if (abs(a) < EPSILON || abs(b) < EPSILON) return 0.0;

    double term1 = -b / (2.0 * (a - b) * (c - b));
    double term2 = -b * b * (-log(abs(a)) + log(abs(b))) / (2.0 * (a - b) * (a - b) * (c - b));

    return term1 + term2;
}

double triangle_integral_1overD(double D1, double D2, double D3, double Area) {
    // Check for zeros
    if (abs(D1) < EPSILON || abs(D2) < EPSILON || abs(D3) < EPSILON) return 0.0;

    // Check if D crosses zero (pole inside triangle)
    double Dmin = min({D1, D2, D3});
    double Dmax = max({D1, D2, D3});
    if (Dmin < -EPSILON && Dmax > EPSILON) return 0.0;

    // All equal case: I = A/D for triangle
    if (abs(D1 - D2) < EPSILON && abs(D2 - D3) < EPSILON) {
        return Area / D1;
    }

    // Two equal case: D1 ≈ D2 ≠ D3
    if (abs(D1 - D2) < EPSILON && abs(D2 - D3) > EPSILON) {
        double a = D1;
        double c = D3;
        if (abs(a - c) < EPSILON) return 0.0;
        double term1 = (a * a - c * c);
        double term2 = -2.0 * a * c * log(abs(a)) + 2.0 * a * c * log(abs(c));
        return (term1 + term2) / (2.0 * (a - c) * (a - c)) * Area;
    }

    // All different - general formula using phi
    double omega1 = phi_function(D1, D2, D3) + phi_function(D1, D3, D2);
    double omega2 = phi_function(D2, D3, D1) + phi_function(D2, D1, D3);
    double omega3 = phi_function(D3, D1, D2) + phi_function(D3, D2, D1);

    return Area * (omega1 + omega2 + omega3);
}

// Compute contribution using simple midpoint rule for now
// This is a fallback - proper implementation would use full tetrahedron splitting
double triangle_contribution(double eps1, double eps2, double eps3,
                             double D1, double D2, double D3,
                             double Area, double eps_F) {
    // Count how many vertices are occupied
    int n_occ = (eps1 < eps_F ? 1 : 0) + (eps2 < eps_F ? 1 : 0) + (eps3 < eps_F ? 1 : 0);

    if (n_occ == 0) return 0.0;  // All unoccupied
    if (n_occ == 3) {
        // All occupied - use full triangle
        return triangle_integral_1overD(D1, D2, D3, Area);
    }

    // Partially occupied - use fraction
    double frac = n_occ / 3.0;
    return triangle_integral_1overD(D1, D2, D3, Area * frac);
}

// Recursive hybrid tetrahedron method for 2D response function
double response_recursive(Vec q, int Nk_coarse, int n_refine) {
    double sum = 0.0;
    double dk = 2.0 * M_PI / Nk_coarse;
    double Area_tri = 0.5 * dk * dk;  // Area of one triangle

    // Loop over coarse grid squares
    for (int i = 0; i < Nk_coarse; i++) {
        for (int j = 0; j < Nk_coarse; j++) {
            // Each square divided into 2 triangles
            Vec k00 = get_k(i, j, Nk_coarse);
            Vec k10 = get_k(i+1, j, Nk_coarse);
            Vec k01 = get_k(i, j+1, Nk_coarse);
            Vec k11 = get_k(i+1, j+1, Nk_coarse);

            // === Triangle 1: (i,j), (i+1,j), (i,j+1) ===
            {
                double eps1 = epsilon(k00);
                double eps2 = epsilon(k10);
                double eps3 = epsilon(k01);

                Vec kq00(k00.x + q.x, k00.y + q.y);
                Vec kq10(k10.x + q.x, k10.y + q.y);
                Vec kq01(k01.x + q.x, k01.y + q.y);

                double eps1q = epsilon(kq00);
                double eps2q = epsilon(kq10);
                double eps3q = epsilon(kq01);

                // Denominator D = eps(k+q) - eps(k)
                double D1 = eps1q - eps1;
                double D2 = eps2q - eps2;
                double D3 = eps3q - eps3;

                // Compute contribution with Fermi surface check
                double contrib = triangle_contribution(eps1, eps2, eps3, D1, D2, D3, Area_tri, mu);
                sum += contrib;
            }

            // === Triangle 2: (i+1,j), (i+1,j+1), (i,j+1) ===
            {
                double eps1 = epsilon(k10);
                double eps2 = epsilon(k11);
                double eps3 = epsilon(k01);

                Vec kq10(k10.x + q.x, k10.y + q.y);
                Vec kq11(k11.x + q.x, k11.y + q.y);
                Vec kq01(k01.x + q.x, k01.y + q.y);

                double eps1q = epsilon(kq10);
                double eps2q = epsilon(kq11);
                double eps3q = epsilon(kq01);

                double D1 = eps1q - eps1;
                double D2 = eps2q - eps2;
                double D3 = eps3q - eps3;

                double contrib = triangle_contribution(eps1, eps2, eps3, D1, D2, D3, Area_tri, mu);
                sum += contrib;
            }
        }
    }

    // Normalize by BZ area (2π)^2
    // Note: Rath & Freeman formula gives positive contribution
    // Response function χ(q) = ∫ [f(k) - f(k+q)] / [ε(k) - ε(k+q)] d²k/(2π)²
    return sum / (4.0 * M_PI * M_PI);
}

// Reference: brute force numerical integration
double response_numerical(Vec q, int Nk) {
    double sum = 0;
    for (int i = 0; i < Nk; i++) {
        for (int j = 0; j < Nk; j++) {
            Vec k = get_k(i, j, Nk);
            double e_k = epsilon(k);
            Vec k_plus_q(k.x + q.x, k.y + q.y);
            double e_qk = epsilon(k_plus_q);
            double f_k = fermi_dirac(e_k);
            double f_qk = fermi_dirac(e_qk);
            double denom = e_k - e_qk;
            if (abs(denom) > 1e-6) {
                sum += -(f_k - f_qk) / denom / (Nk * Nk);
            }
        }
    }
    return sum;
}

int main() {
    cout << "Testing Recursive Hybrid Tetrahedron Method for χ(q)\n";
    cout << "====================================================\n\n";
    cout << "Model: 2D tight-binding, t=" << t << ", mu=" << mu << "\n\n";

    // Test on 3 q-points as requested
    vector<Vec> q_points = {
        Vec(0.2, 0.0),
        Vec(0.5, 0.5),
        Vec(M_PI/4, M_PI/4)
    };

    // Reference: high-resolution numerical
    int Nk_ref = 300;
    cout << "Computing reference values (Nk = " << Nk_ref << ")...\n";
    vector<double> ref_vals;
    for (const auto& q : q_points) {
        double ref = response_numerical(q, Nk_ref);
        ref_vals.push_back(ref);
        cout << "  q=(" << fixed << setprecision(3) << q.x << ", " << q.y << "): "
             << setprecision(6) << ref << "\n";
    }

    // Test recursive method with different grid sizes
    vector<int> Nk_test = {10, 20, 30, 50, 60, 80, 300};

    cout << "\n\nRecursive Tetrahedron Method Results:\n";
    cout << "--------------------------------------\n";
    cout << fixed << setprecision(6);

    for (size_t iq = 0; iq < q_points.size(); iq++) {
        cout << "\nq = (" << q_points[iq].x << ", " << q_points[iq].y << ")\n";
        cout << "Reference: " << setw(12) << ref_vals[iq] << "\n";
        cout << "Nk   χ(q)         Error(%)\n";

        for (int Nk : Nk_test) {
            double chi_rec = response_recursive(q_points[iq], Nk, 0);
            double error = abs(chi_rec - ref_vals[iq]) / (abs(ref_vals[iq]) + 1e-10) * 100.0;
            cout << setw(3) << Nk << "  " << setw(12) << chi_rec
                 << "  " << setw(7) << error << "\n";
        }
    }

    cout << "\n\nTarget: Achieve ~20% error margin\n";

    return 0;
}
