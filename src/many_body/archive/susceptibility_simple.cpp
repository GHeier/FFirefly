/*
 * Simplified implementation comparing tetrahedron-inspired integration
 * with standard numerical integration for generalized susceptibility χ(q)
 *
 * Reference: J. Rath and A. J. Freeman, Phys. Rev. B 11, 2109 (1975)
 *
 * This version implements the core analytical formulas (Eq. 17-21) for
 * integration over a tetrahedron with linearized energy denominator.
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <fstream>

const double PI = 3.14159265358979323846;
const double EPSILON = 1e-10;

// 3D vector
struct Vec3 {
    double x, y, z;
    Vec3(double x_=0, double y_=0, double z_=0) : x(x_), y(y_), z(z_) {}
    Vec3 operator+(const Vec3& v) const { return Vec3(x+v.x, y+v.y, z+v.z); }
    Vec3 operator-(const Vec3& v) const { return Vec3(x-v.x, y-v.y, z-v.z); }
    Vec3 operator*(double s) const { return Vec3(x*s, y*s, z*s); }
    double dot(const Vec3& v) const { return x*v.x + y*v.y + z*v.z; }
};

/*
 * 2D Tight-binding model: E(k) = -2t[cos(kx*a) + cos(ky*a)]
 */
class TightBinding2D {
public:
    double t, a;
    TightBinding2D(double t_=1.0, double a_=1.0) : t(t_), a(a_) {}

    double energy(double kx, double ky) const {
        return -2.0 * t * (std::cos(kx*a) + std::cos(ky*a));
    }

    double fermi_energy() const { return 0.0; } // half-filling
};

/*
 * Core integral formula for triangle (2D) with linearized denominator
 * For a triangle with values V0, V1, V2 at corners and area A:
 * I = ∫_triangle dA / V(r)
 *
 * This is adapted from the tetrahedron formulas
 */
double triangle_integral(double V0, double V1, double V2, double Area) {
    // Sort values
    std::vector<double> V = {V0, V1, V2};
    std::sort(V.begin(), V.end(), std::greater<double>());
    // V[0] ≥ V[1] ≥ V[2]

    // Check all have same sign
    if ((V[0] > EPSILON && V[2] < -EPSILON) || (V[0] < -EPSILON && V[2] > EPSILON)) {
        // Values cross zero - integral is problematic
        return 0.0;
    }

    // Check for zeros
    if (std::abs(V[0]) < EPSILON || std::abs(V[1]) < EPSILON || std::abs(V[2]) < EPSILON) {
        return 0.0;
    }

    // All equal
    if (std::abs(V[0]-V[1]) < EPSILON && std::abs(V[1]-V[2]) < EPSILON) {
        return Area / V[0];
    }

    // V[0] = V[1] != V[2]
    if (std::abs(V[0]-V[1]) < EPSILON && std::abs(V[1]-V[2]) > EPSILON) {
        double Va = V[0], Vb = V[2];
        double dV = Va - Vb;
        if (std::abs(dV) < EPSILON) return 0.0;
        return 2.0 * Area * (Va/(dV*dV) * std::log(std::abs(Va/Vb)) + 1.0/dV);
    }

    // All different - use simplified formula for triangle
    double denom1 = (V[0]-V[2])*(V[0]-V[1]);
    double denom2 = (V[1]-V[2])*(V[1]-V[0]);

    if (std::abs(denom1*denom2) < EPSILON) return 0.0;

    double term1 = V[0] / denom1 * std::log(std::abs(V[0]/V[2]));
    double term2 = V[1] / denom2 * std::log(std::abs(V[1]/V[2]));

    return 2.0 * Area * (term1 + term2);
}

/*
 * Core tetrahedron integral formula from Rath & Freeman Eq. (17)
 * (Kept for reference but we'll use triangle_integral for 2D)
 *
 * Integral: I = ∫_tetrahedron d³k / V(k)
 * where V(k) is linearized as V(k) = V4 + coeffs·(k-k4)
 *
 * Returns: 3Ω * [V1²/D1 * ln|V1/V4| + V2²/D2 * ln|V2/V4| + V3²/D3 * ln|V3/V4|]
 * where Ω = volume of tetrahedron
 *       V1 ≥ V2 ≥ V3 ≥ V4 are values at corners
 *       Di = (Vi-V4)(Vi-V3)(Vi-V2) for i=1,2,3
 */
double tetrahedron_integral(double V1, double V2, double V3, double V4, double Omega) {
    // Handle degenerate cases first (Eq. 18-21)

    // Check for zeros
    if (std::abs(V1) < EPSILON && std::abs(V2) < EPSILON &&
        std::abs(V3) < EPSILON && std::abs(V4) < EPSILON) {
        return 0.0; // All zeros - singular
    }

    // Eq. (18): All equal V1=V2=V3=V4=V ≠ 0
    if (std::abs(V1-V2) < EPSILON && std::abs(V2-V3) < EPSILON &&
        std::abs(V3-V4) < EPSILON) {
        return Omega / V1;
    }

    // Eq. (19): V1=V2=V3=V ≠ V4
    if (std::abs(V1-V2) < EPSILON && std::abs(V2-V3) < EPSILON &&
        std::abs(V3-V4) > EPSILON) {
        double V = V1;
        if (std::abs(V-V4) < EPSILON) return 0.0;
        double term1 = V4*V4 / std::pow(V-V4, 3) * std::log(std::abs(V/V4));
        double term2 = (0.5*V*V + 0.25*V4*V4 - V*V4) / std::pow(V-V4, 3);
        return 3.0 * Omega * (term1 + term2);
    }

    // Eq. (20): V1=V2=V ≠ V3=V4=V'
    if (std::abs(V1-V2) < EPSILON && std::abs(V3-V4) < EPSILON &&
        std::abs(V1-V3) > EPSILON) {
        double V = V1, Vp = V3;
        if (std::abs(V-Vp) < EPSILON) return 0.0;
        double term1 = V*Vp / std::pow(V-Vp, 3) * std::log(std::abs(V/Vp));
        double term2 = (V + Vp) / std::pow(V-Vp, 2);
        return 3.0 * Omega * (term1 + term2);
    }

    // Eq. (21): V1=V2 ≠ V3 ≠ V4
    if (std::abs(V1-V2) < EPSILON && std::abs(V1-V3) > EPSILON &&
        std::abs(V3-V4) > EPSILON) {
        double term1 = V2*V2 / ((V2-V4)*(V2-V4)*(V2-V3)*(V2-V3)) * std::log(std::abs(V2/V4));
        double term2 = V3*V3 / ((V4-V3)*(V4-V3)*(V2-V3)*(V2-V3)) * std::log(std::abs(V3/V4));
        double term3 = V3 / ((V2-V3)*(V4-V3));
        return 3.0 * Omega * (term1 + term2 + term3);
    }

    // Eq. (17): General case - all different
    double D1 = (V1-V4)*(V1-V3)*(V1-V2);
    double D2 = (V2-V4)*(V2-V3)*(V2-V1);
    double D3 = (V3-V4)*(V3-V2)*(V3-V1);

    if (std::abs(D1*D2*D3) < EPSILON) {
        return 0.0; // Singular
    }

    // Check if V crosses zero - need special handling
    bool has_positive = (V1 > EPSILON) || (V2 > EPSILON) || (V3 > EPSILON) || (V4 > EPSILON);
    bool has_negative = (V1 < -EPSILON) || (V2 < -EPSILON) || (V3 < -EPSILON) || (V4 < -EPSILON);

    if (has_positive && has_negative) {
        // V changes sign across tetrahedron
        // The integral has a singularity and needs careful treatment
        // For now, return 0 to avoid NaN (proper treatment requires subdivision)
        return 0.0;
    }

    if (std::abs(V1) < EPSILON || std::abs(V2) < EPSILON ||
        std::abs(V3) < EPSILON || std::abs(V4) < EPSILON) {
        // Has value very close to zero - risky for log
        return 0.0;
    }

    // Check that all values have same sign for safe logarithm
    if ((V1*V4 < 0) || (V2*V4 < 0) || (V3*V4 < 0)) {
        return 0.0; // Different signs - avoid log of negative number
    }

    double term1 = V1*V1 / D1 * std::log(std::abs(V1/V4));
    double term2 = V2*V2 / D2 * std::log(std::abs(V2/V4));
    double term3 = V3*V3 / D3 * std::log(std::abs(V3/V4));

    return 3.0 * Omega * (term1 + term2 + term3);
}

/*
 * Calculate χ(q) using tetrahedron-based method
 * We divide BZ into rectangles, each split into 2 triangles
 */
double chi_tetrahedron_method(const Vec3& q, const TightBinding2D& tb, int Nk) {
    double Ef = tb.fermi_energy();
    double chi = 0.0;

    double dk = 2.0*PI / (tb.a * Nk);
    double triangle_area = 0.5 * dk * dk; // Area of each triangle in 2D

    int contributions = 0;

    for (int ix = 0; ix < Nk; ix++) {
        for (int iy = 0; iy < Nk; iy++) {
            double kx0 = -PI/tb.a + ix * dk;
            double ky0 = -PI/tb.a + iy * dk;

            // Each square divided into 2 triangles (tetrahedrons in 2D)
            // Triangle 1: (kx0,ky0), (kx0+dk,ky0), (kx0,ky0+dk)
            {
                Vec3 k0(kx0, ky0, 0);
                Vec3 k1(kx0+dk, ky0, 0);
                Vec3 k2(kx0, ky0+dk, 0);

                // Energies at corners
                double E0 = tb.energy(k0.x, k0.y);
                double E1 = tb.energy(k1.x, k1.y);
                double E2 = tb.energy(k2.x, k2.y);

                // Energies at k+q
                double Ep0 = tb.energy(k0.x+q.x, k0.y+q.y);
                double Ep1 = tb.energy(k1.x+q.x, k1.y+q.y);
                double Ep2 = tb.energy(k2.x+q.x, k2.y+q.y);

                // Check if triangle spans Fermi surface appropriately
                double Emin = std::min({E0, E1, E2});
                double Emax = std::max({E0, E1, E2});
                double Epmin = std::min({Ep0, Ep1, Ep2});
                double Epmax = std::max({Ep0, Ep1, Ep2});

                // Need E < Ef < E' for contribution
                if (Emin < Ef && Epmax > Ef) {
                    // Calculate V_i = E'(k_i+q) - E(k_i)
                    double V0 = Ep0 - E0;
                    double V1 = Ep1 - E1;
                    double V2 = Ep2 - E2;

                    double integral = triangle_integral(V0, V1, V2, triangle_area);
                    chi += integral;
                    if (std::abs(integral) > EPSILON) contributions++;
                }
            }

            // Triangle 2: (kx0+dk,ky0), (kx0+dk,ky0+dk), (kx0,ky0+dk)
            {
                Vec3 k0(kx0+dk, ky0, 0);
                Vec3 k1(kx0+dk, ky0+dk, 0);
                Vec3 k2(kx0, ky0+dk, 0);

                double E0 = tb.energy(k0.x, k0.y);
                double E1 = tb.energy(k1.x, k1.y);
                double E2 = tb.energy(k2.x, k2.y);

                double Ep0 = tb.energy(k0.x+q.x, k0.y+q.y);
                double Ep1 = tb.energy(k1.x+q.x, k1.y+q.y);
                double Ep2 = tb.energy(k2.x+q.x, k2.y+q.y);

                double Emin = std::min({E0, E1, E2});
                double Emax = std::max({E0, E1, E2});
                double Epmin = std::min({Ep0, Ep1, Ep2});
                double Epmax = std::max({Ep0, Ep1, Ep2});

                if (Emin < Ef && Epmax > Ef) {
                    double V0 = Ep0 - E0;
                    double V1 = Ep1 - E1;
                    double V2 = Ep2 - E2;

                    double integral = triangle_integral(V0, V1, V2, triangle_area);
                    chi += integral;
                    if (std::abs(integral) > EPSILON) contributions++;
                }
            }
        }
    }

    // Multiply by spin factor
    chi *= 2.0;

    return chi;
}

/*
 * Standard numerical integration
 */
double chi_numerical(const Vec3& q, const TightBinding2D& tb, int Nk) {
    double Ef = tb.fermi_energy();
    double chi = 0.0;
    double dk = 2.0*PI / (tb.a * Nk);

    for (int ix = 0; ix < Nk; ix++) {
        for (int iy = 0; iy < Nk; iy++) {
            double kx = -PI/tb.a + (ix + 0.5) * dk;
            double ky = -PI/tb.a + (iy + 0.5) * dk;

            double Ek = tb.energy(kx, ky);
            double Ekq = tb.energy(kx+q.x, ky+q.y);

            // f(Ek)[1-f(Ekq)] at T=0
            if (Ek < Ef && Ekq > Ef) {
                double denom = Ekq - Ek;
                if (std::abs(denom) > EPSILON) {
                    chi += 1.0 / denom;
                }
            }
        }
    }

    // dk*dk for area element, factor of 2 for spin
    chi *= dk * dk * 2.0;

    return chi;
}

int main() {
    std::cout << "==============================================\n";
    std::cout << "Generalized Susceptibility χ(q) Calculation\n";
    std::cout << "Tetrahedron Method vs Numerical Integration\n";
    std::cout << "==============================================\n\n";

    TightBinding2D tb(1.0, 1.0);
    std::cout << "2D Tight-Binding Model: E(k) = -2t[cos(kx) + cos(ky)]\n";
    std::cout << "  t = " << tb.t << ", a = " << tb.a << "\n";
    std::cout << "  Ef = " << tb.fermi_energy() << " (half-filling)\n\n";

    int Nk = 60;
    std::cout << "k-mesh: " << Nk << " × " << Nk << "\n\n";

    std::cout << std::setw(12) << "qx/π"
              << std::setw(16) << "χ_tetrahedron"
              << std::setw(16) << "χ_numerical"
              << std::setw(14) << "ratio\n";
    std::cout << std::string(58, '-') << "\n";

    std::ofstream outfile("chi_comparison.dat");
    outfile << "# qx/pi  chi_tetrahedron  chi_numerical  ratio\n";

    for (int iq = 0; iq <= 20; iq++) {
        double qx = (iq * PI / tb.a) / 20.0;
        Vec3 q(qx, 0, 0);

        double chi_tet = chi_tetrahedron_method(q, tb, Nk);
        double chi_num = chi_numerical(q, tb, Nk);
        double ratio = (std::abs(chi_num) > EPSILON) ? chi_tet/chi_num : 0.0;

        std::cout << std::setw(12) << std::fixed << std::setprecision(4) << qx/PI
                  << std::setw(16) << std::setprecision(6) << chi_tet
                  << std::setw(16) << chi_num
                  << std::setw(14) << std::setprecision(4) << ratio << "\n";

        outfile << qx/PI << "  " << chi_tet << "  " << chi_num << "  " << ratio << "\n";
    }

    outfile.close();
    std::cout << "\nResults saved to chi_comparison.dat\n";

    return 0;
}
