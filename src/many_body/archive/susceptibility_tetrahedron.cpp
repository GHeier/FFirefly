/*
 * Implementation of the analytic tetrahedron linear energy method
 * for calculating generalized magnetic susceptibility χ(q)
 *
 * Reference: J. Rath and A. J. Freeman, Phys. Rev. B 11, 2109 (1975)
 * "Generalized magnetic susceptibilities in metals: Application of the
 *  analytic tetrahedron linear energy method to Sc"
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <algorithm>
#include <iomanip>
#include <fstream>

const double PI = 3.14159265358979323846;
const double EPSILON = 1e-10;

// 3D vector class
struct Vec3 {
    double x, y, z;

    Vec3() : x(0), y(0), z(0) {}
    Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    Vec3 operator+(const Vec3& v) const { return Vec3(x+v.x, y+v.y, z+v.z); }
    Vec3 operator-(const Vec3& v) const { return Vec3(x-v.x, y-v.y, z-v.z); }
    Vec3 operator*(double s) const { return Vec3(x*s, y*s, z*s); }
    double dot(const Vec3& v) const { return x*v.x + y*v.y + z*v.z; }
    Vec3 cross(const Vec3& v) const {
        return Vec3(y*v.z - z*v.y, z*v.x - x*v.z, x*v.y - y*v.x);
    }
    double norm() const { return std::sqrt(x*x + y*y + z*z); }
};

// Tetrahedron with 4 corners
struct Tetrahedron {
    Vec3 k[4];  // k-points at corners
    double E[4]; // energies at corners (band n)
    double Ep[4]; // energies at corners (band n' at k+q)

    double volume() const {
        Vec3 v1 = k[1] - k[0];
        Vec3 v2 = k[2] - k[0];
        Vec3 v3 = k[3] - k[0];
        return std::abs(v1.dot(v2.cross(v3))) / 6.0;
    }
};

// Sort energies and return permutation indices
std::vector<int> argsort(const double* arr, int n) {
    std::vector<int> idx(n);
    for (int i = 0; i < n; i++) idx[i] = i;
    std::sort(idx.begin(), idx.end(), [&arr](int i, int j) { return arr[i] < arr[j]; });
    return idx;
}

/*
 * Calculate contribution to χ(q) from a single tetrahedron
 * Using equations (17)-(21) from Rath & Freeman (1975)
 *
 * The integral is:
 * I = ∫ d³k / [E_n'(k+q) - E_n(k)]
 *
 * With linearized energy denominator V_i = E_n'(k_i+q) - E_n(k_i)
 */
double tetrahedron_chi_contribution(const Tetrahedron& tet, double Ef) {
    // For constant matrix element approximation, we need to check
    // if the tetrahedron has any region where f(E_n)[1-f(E_n')] ≠ 0
    // This requires geometric analysis of Fermi surface intersections.

    // For simplicity in this implementation, we use a constant-matrix-element
    // approximation and integrate over the full tetrahedron, then weight
    // by occupation. A more sophisticated implementation would subdivide
    // the tetrahedron based on Fermi surface intersections (see Figs 1-3 in paper).

    // Check if tetrahedron has contributing states
    // We need f(E)[1-f(E')] ≠ 0, i.e., SOME occupied states (E < Ef)
    // AND SOME empty shifted states (E' > Ef)

    // Count occupied and empty states at corners
    int n_occupied = 0, n_empty_shifted = 0;
    static int debug_count = 0;
    for (int i = 0; i < 4; i++) {
        if (tet.E[i] < Ef) n_occupied++;
        if (tet.Ep[i] > Ef) n_empty_shifted++;
    }

    if (debug_count < 5 && (n_occupied > 0 || n_empty_shifted < 4)) {
        std::cout << "DEBUG contrib: n_occ=" << n_occupied << " n_empty=" << n_empty_shifted
                  << " E=[" << tet.E[0] << "," << tet.E[1] << "," << tet.E[2] << "," << tet.E[3] << "]"
                  << " Ep=[" << tet.Ep[0] << "," << tet.Ep[1] << "," << tet.Ep[2] << "," << tet.Ep[3] << "]\n";
        debug_count++;
    }

    // Need at least one occupied AND one empty shifted state
    if (n_occupied == 0 || n_empty_shifted == 0) return 0.0;

    // Calculate V_i = E_n'(k_i+q) - E_n(k_i)
    double V[4];
    for (int i = 0; i < 4; i++) {
        V[i] = tet.Ep[i] - tet.E[i];
    }

    // Check if all V values have the same sign - if so, no interband transitions
    bool all_positive = true, all_negative = true;
    for (int i = 0; i < 4; i++) {
        if (V[i] < 0) all_positive = false;
        if (V[i] > 0) all_negative = false;
    }
    if (all_positive || all_negative) {
        // For χ(q), we need transitions between occupied and empty states
        // If V doesn't change sign, integral contribution may be suppressed
        // but we'll still compute it
    }

    // Sort V values: V4 ≤ V3 ≤ V2 ≤ V1
    std::vector<int> idx = argsort(V, 4);
    double Vs[4];
    for (int i = 0; i < 4; i++) {
        Vs[i] = V[idx[3-i]]; // Reverse order for decreasing sort
    }
    double V1 = Vs[0], V2 = Vs[1], V3 = Vs[2], V4 = Vs[3];

    double Omega = tet.volume();

    // For 2D (degenerate tetrahedrons), volume formula gives 0
    // Use triangle area instead: 0.5 * |edge1 × edge2|
    if (Omega < EPSILON) {
        Vec3 edge1 = tet.k[1] - tet.k[0];
        Vec3 edge2 = tet.k[2] - tet.k[0];
        Vec3 cross_prod = Vec3(
            edge1.y * edge2.z - edge1.z * edge2.y,
            edge1.z * edge2.x - edge1.x * edge2.z,
            edge1.x * edge2.y - edge1.y * edge2.x
        );
        Omega = 0.5 * std::sqrt(cross_prod.x*cross_prod.x +
                                cross_prod.y*cross_prod.y +
                                cross_prod.z*cross_prod.z);
    }

    if (Omega < EPSILON) return 0.0;

    // Apply equations (17)-(21) for different cases
    double result = 0.0;

    // Case (i): V1 = V2 = V3 = V4 = V ≠ 0
    if (std::abs(V1 - V4) < EPSILON && std::abs(V1) > EPSILON) {
        // Eq. (18)
        result = Omega / V1;
    }
    // Case (ii): V1 = V2 = V3 = V ≠ V4
    else if (std::abs(V1 - V3) < EPSILON && std::abs(V1 - V4) > EPSILON) {
        // Eq. (19)
        double V = V1;
        if (std::abs(V - V4) > EPSILON) {
            result = 3.0 * Omega * (
                V4*V4 / ((V - V4) * (V - V4) * (V - V4)) * std::log(std::abs(V/V4))
                + 0.5 * (V*V + 0.5*V4*V4 - 2.0*V*V4) / ((V - V4) * (V - V4) * (V - V4))
            );
        }
    }
    // Case (iii): V1 = V2 = V ≠ V3 = V4 = V'
    else if (std::abs(V1 - V2) < EPSILON && std::abs(V3 - V4) < EPSILON &&
             std::abs(V1 - V3) > EPSILON) {
        // Eq. (20)
        double V = V1, Vp = V3;
        if (std::abs(V - Vp) > EPSILON) {
            result = 3.0 * Omega * (
                V*Vp / ((V - Vp) * (V - Vp) * (V - Vp)) * std::log(std::abs(V/Vp))
                + (V + Vp) / ((V - Vp) * (V - Vp))
            );
        }
    }
    // Case (iv): V1 = V2 ≠ V ≠ V3 ≠ V4 (general case)
    else if (std::abs(V1 - V2) < EPSILON && std::abs(V1 - V3) > EPSILON &&
             std::abs(V1 - V4) > EPSILON && std::abs(V3 - V4) > EPSILON) {
        // Eq. (21)
        result = 3.0 * Omega * (
            V2*V2 / ((V2 - V4) * (V2 - V4) * (V2 - V3) * (V2 - V3)) * std::log(std::abs(V2/V4))
            + V3*V3 / ((V4 - V3) * (V4 - V3) * (V2 - V3) * (V2 - V3)) * std::log(std::abs(V3/V4))
            + V3 / ((V2 - V3) * (V4 - V3))
        );
    }
    // Case: All different (Eq. 17 - general formula)
    else {
        // Check if all V values are distinct
        bool all_different = true;
        for (int i = 0; i < 3; i++) {
            for (int j = i+1; j < 4; j++) {
                if (std::abs(Vs[i] - Vs[j]) < EPSILON) {
                    all_different = false;
                    break;
                }
            }
            if (!all_different) break;
        }

        if (all_different && std::abs(V1*V2*V3*V4) > EPSILON) {
            // Eq. (17)
            double D1 = (V1 - V4) * (V1 - V3) * (V1 - V2);
            double D2 = (V2 - V4) * (V2 - V3) * (V2 - V1);
            double D3 = (V3 - V4) * (V3 - V2) * (V3 - V1);

            if (std::abs(D1*D2*D3) > EPSILON) {
                result = 3.0 * Omega * (
                    V1*V1 / D1 * std::log(std::abs(V1/V4))
                    + V2*V2 / D2 * std::log(std::abs(V2/V4))
                    + V3*V3 / D3 * std::log(std::abs(V3/V4))
                );
            }
        }
    }

    return result;
}

/*
 * 2D Tight-binding model on square lattice
 * E(k) = -2t[cos(k_x) + cos(k_y)] - mu
 */
class TightBinding2D {
public:
    double t;  // hopping parameter
    double a;  // lattice constant
    double mu; // chemical potential

    TightBinding2D(double t_ = 1.0, double a_ = 1.0, double mu_ = 0.0)
        : t(t_), a(a_), mu(mu_) {}

    // Energy dispersion
    double energy(const Vec3& k) const {
        return -2.0 * t * (std::cos(k.x * a) + std::cos(k.y * a)) - mu;
    }

    // Fermi energy (always at 0 with chemical potential included in energy)
    double fermi_energy() const {
        return 0.0;
    }
};

/*
 * Create mesh of tetrahedrons for 2D square lattice
 * Divide each square into 2 triangles (tetrahedrons in 2D are triangles)
 */
std::vector<Tetrahedron> create_2D_mesh(int Nk, double Lx, double Ly) {
    std::vector<Tetrahedron> tetrahedra;

    double dx = Lx / Nk;
    double dy = Ly / Nk;

    // For 2D, we use degenerate tetrahedrons (all kz = 0)
    // Each square is divided into 2 triangles
    for (int ix = 0; ix < Nk; ix++) {
        for (int iy = 0; iy < Nk; iy++) {
            double x0 = ix * dx - Lx/2;
            double y0 = iy * dy - Ly/2;

            // Triangle 1: (x0,y0), (x0+dx,y0), (x0,y0+dy)
            Tetrahedron tet1;
            tet1.k[0] = Vec3(x0, y0, 0);
            tet1.k[1] = Vec3(x0 + dx, y0, 0);
            tet1.k[2] = Vec3(x0, y0 + dy, 0);
            tet1.k[3] = Vec3(x0, y0, 0); // Duplicate for 4th corner (degenerate)
            tetrahedra.push_back(tet1);

            // Triangle 2: (x0+dx,y0), (x0+dx,y0+dy), (x0,y0+dy)
            Tetrahedron tet2;
            tet2.k[0] = Vec3(x0 + dx, y0, 0);
            tet2.k[1] = Vec3(x0 + dx, y0 + dy, 0);
            tet2.k[2] = Vec3(x0, y0 + dy, 0);
            tet2.k[3] = Vec3(x0 + dx, y0, 0); // Duplicate
            tetrahedra.push_back(tet2);
        }
    }

    return tetrahedra;
}

/*
 * Calculate χ(q) using tetrahedron method
 */
double chi_tetrahedron(const Vec3& q, const TightBinding2D& tb,
                       const std::vector<Tetrahedron>& mesh) {
    double Ef = tb.fermi_energy();
    double chi = 0.0;

    int contributions = 0;
    bool debug_first = true;

    for (auto tet : mesh) {
        // Calculate energies at corners
        for (int i = 0; i < 4; i++) {
            tet.E[i] = tb.energy(tet.k[i]);
            Vec3 kpq = tet.k[i] + q;
            tet.Ep[i] = tb.energy(kpq);
        }

        if (debug_first) {
            std::cout << "DEBUG first tetrahedron:\n";
            for (int i = 0; i < 4; i++) {
                std::cout << "  k[" << i << "] = (" << tet.k[i].x << ", " << tet.k[i].y
                          << ") E=" << tet.E[i] << " Ep=" << tet.Ep[i] << "\n";
            }
            debug_first = false;
        }

        double contrib = tetrahedron_chi_contribution(tet, Ef);
        chi += contrib;
        if (std::abs(contrib) > EPSILON) contributions++;
    }

    std::cout << "  Contributions: " << contributions << " / " << mesh.size() << "\n";

    // Normalize by BZ volume (factor of 2 for spin)
    double BZ_area = 4.0 * PI * PI / (tb.a * tb.a);
    chi *= 2.0 / BZ_area;

    return chi;
}

/*
 * Calculate χ(q) using simple numerical integration (for comparison)
 */
double chi_numerical(const Vec3& q, const TightBinding2D& tb, int Nk) {
    double Ef = tb.fermi_energy();
    double chi = 0.0;

    double dk = 2.0 * PI / (tb.a * Nk);
    int count = 0;

    for (int ix = 0; ix < Nk; ix++) {
        for (int iy = 0; iy < Nk; iy++) {
            double kx = -PI/tb.a + (ix + 0.5) * dk;
            double ky = -PI/tb.a + (iy + 0.5) * dk;

            Vec3 k(kx, ky, 0);
            Vec3 kpq = k + q;

            double Ek = tb.energy(k);
            double Ekq = tb.energy(kpq);

            // Check Fermi factors: f(Ek)[1-f(Ekq)]
            bool occupied = (Ek < Ef);
            bool unoccupied = (Ekq > Ef);

            if (occupied && unoccupied) {
                double denom = Ekq - Ek;
                if (std::abs(denom) > EPSILON) {
                    chi += 1.0 / denom;
                    count++;
                }
            }
        }
    }

    // Multiply by d³k volume element and spin factor
    double dV = dk * dk;
    chi *= dV * 2.0;

    return chi;
}

int main() {
    std::cout << "========================================\n";
    std::cout << "Susceptibility Calculation\n";
    std::cout << "Tetrahedron Method vs Numerical Integration\n";
    std::cout << "========================================\n\n";

    // Setup 2D tight-binding model
    double t = 1.0;
    double a = 1.0;
    double mu = -1.0;  // Chemical potential
    TightBinding2D tb(t, a, mu);

    std::cout << "2D Tight-Binding Model Parameters:\n";
    std::cout << "  Hopping t = " << t << "\n";
    std::cout << "  Lattice constant a = " << a << "\n";
    std::cout << "  Chemical potential mu = " << mu << "\n";
    std::cout << "  Fermi energy Ef = " << tb.fermi_energy() << "\n\n";

    // Create k-mesh
    int Nk = 40;
    double Lx = 2.0 * PI / a;
    double Ly = 2.0 * PI / a;

    std::cout << "Creating k-mesh with Nk = " << Nk << "...\n";
    auto mesh = create_2D_mesh(Nk, Lx, Ly);
    std::cout << "Created " << mesh.size() << " tetrahedra\n\n";

    // Calculate χ(q) along (π,π) direction
    std::cout << "Calculating χ(q) along (π,π) direction...\n\n";
    std::cout << std::setw(10) << "q/(π,π)"
              << std::setw(15) << "χ_tetra"
              << std::setw(15) << "χ_numerical"
              << std::setw(15) << "difference\n";
    std::cout << std::string(55, '-') << "\n";

    std::ofstream outfile("susceptibility_comparison.dat");
    outfile << "# q/(pi,pi)  chi_tetrahedron  chi_numerical  difference\n";

    int Nq = 20;
    for (int iq = 0; iq <= Nq; iq++) {
        double q_mag = (iq * PI / a) / Nq;
        Vec3 q(q_mag, q_mag, 0);  // Along (π,π) direction

        double chi_tet = chi_tetrahedron(q, tb, mesh);
        double chi_num = chi_numerical(q, tb, Nk);
        double diff = std::abs(chi_tet - chi_num);

        std::cout << std::setw(10) << std::fixed << std::setprecision(4) << q_mag/PI
                  << std::setw(15) << std::setprecision(6) << chi_tet
                  << std::setw(15) << chi_num
                  << std::setw(15) << diff << "\n";

        outfile << q_mag/PI << "  " << chi_tet << "  " << chi_num << "  " << diff << "\n";
        std::cout << q_mag/PI << "  " << chi_tet << "  " << chi_num << "  " << diff << std::endl;
    }

    outfile.close();
    std::cout << "\nResults written to susceptibility_comparison.dat\n";

    return 0;
}
