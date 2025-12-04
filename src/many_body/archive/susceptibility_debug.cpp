/*
 * DEBUG version with detailed output
 * Testing susceptibility calculation with mu = -1.0 along (pi,pi) direction
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <fstream>

const double PI = 3.14159265358979323846;
const double EPSILON = 1e-10;

struct Vec3 {
    double x, y, z;
    Vec3(double x_=0, double y_=0, double z_=0) : x(x_), y(y_), z(z_) {}
    Vec3 operator+(const Vec3& v) const { return Vec3(x+v.x, y+v.y, z+v.z); }
    Vec3 operator-(const Vec3& v) const { return Vec3(x-v.x, y-v.y, z-v.z); }
};

class TightBinding2D {
public:
    double t, a, mu;
    TightBinding2D(double t_=1.0, double a_=1.0, double mu_=0.0) : t(t_), a(a_), mu(mu_) {}

    double energy(double kx, double ky) const {
        return -2.0 * t * (std::cos(kx*a) + std::cos(ky*a)) - mu;
    }

    double fermi_energy() const { return 0.0; } // Always at 0 with chemical potential in energy
};

// Triangle integral - simplified and debugged version
double triangle_integral(double V0, double V1, double V2, double Area, bool debug=false) {
    std::vector<double> V = {V0, V1, V2};
    std::sort(V.begin(), V.end(), std::greater<double>());
    // V[0] >= V[1] >= V[2]

    if (debug) {
        std::cout << "    Triangle: V = [" << V[0] << ", " << V[1] << ", " << V[2] << "], Area=" << Area << "\n";
    }

    // Check if values cross zero - this means a pole in the integral
    if ((V[0] > EPSILON && V[2] < -EPSILON) || (V[0] < -EPSILON && V[2] > EPSILON)) {
        if (debug) std::cout << "    -> Crosses zero, skipping\n";
        return 0.0;
    }

    // Check for zeros
    if (std::abs(V[0]) < EPSILON || std::abs(V[1]) < EPSILON || std::abs(V[2]) < EPSILON) {
        if (debug) std::cout << "    -> Contains zero, skipping\n";
        return 0.0;
    }

    // All equal
    if (std::abs(V[0]-V[1]) < EPSILON && std::abs(V[1]-V[2]) < EPSILON) {
        double result = Area / V[0];
        if (debug) std::cout << "    -> All equal: I = " << result << "\n";
        return result;
    }

    // V[0] = V[1] != V[2]
    if (std::abs(V[0]-V[1]) < EPSILON && std::abs(V[1]-V[2]) > EPSILON) {
        double Va = V[0], Vb = V[2];
        double dV = Va - Vb;
        if (std::abs(dV) < EPSILON) return 0.0;
        double result = 2.0 * Area * (Va/(dV*dV) * std::log(std::abs(Va/Vb)) + 1.0/dV);
        if (debug) std::cout << "    -> Two equal: I = " << result << "\n";
        return result;
    }

    // All different
    double denom1 = (V[0]-V[2])*(V[0]-V[1]);
    double denom2 = (V[1]-V[2])*(V[1]-V[0]);

    if (std::abs(denom1*denom2) < EPSILON) return 0.0;

    double term1 = V[0] / denom1 * std::log(std::abs(V[0]/V[2]));
    double term2 = V[1] / denom2 * std::log(std::abs(V[1]/V[2]));
    double result = 2.0 * Area * (term1 + term2);

    if (debug) std::cout << "    -> All different: I = " << result << "\n";
    return result;
}

double chi_tetrahedron_method(const Vec3& q, const TightBinding2D& tb, int Nk, bool verbose=false) {
    double Ef = tb.fermi_energy();
    double chi = 0.0;

    double dk = 2.0*PI / (tb.a * Nk);
    double triangle_area = 0.5 * dk * dk;

    int total_triangles = 0;
    int contributing_triangles = 0;

    for (int ix = 0; ix < Nk; ix++) {
        for (int iy = 0; iy < Nk; iy++) {
            double kx0 = -PI/tb.a + ix * dk;
            double ky0 = -PI/tb.a + iy * dk;

            // Triangle 1
            {
                total_triangles++;
                Vec3 k0(kx0, ky0, 0);
                Vec3 k1(kx0+dk, ky0, 0);
                Vec3 k2(kx0, ky0+dk, 0);

                double E0 = tb.energy(k0.x, k0.y);
                double E1 = tb.energy(k1.x, k1.y);
                double E2 = tb.energy(k2.x, k2.y);

                double Ep0 = tb.energy(k0.x+q.x, k0.y+q.y);
                double Ep1 = tb.energy(k1.x+q.x, k1.y+q.y);
                double Ep2 = tb.energy(k2.x+q.x, k2.y+q.y);

                // At T=0: need f(E)[1-f(E')] = 1, i.e., E < Ef AND E' > Ef
                // Check if ANY corner satisfies this
                bool has_contribution = false;
                for (int i = 0; i < 3; i++) {
                    double E_corner = (i==0) ? E0 : (i==1 ? E1 : E2);
                    double Ep_corner = (i==0) ? Ep0 : (i==1 ? Ep1 : Ep2);

                    if (E_corner < Ef || Ep_corner > Ef) {
                        has_contribution = true;
                        break;
                    }
                }

                // Better criterion: check if the integrand has support
                // We need SOME occupied states AND SOME empty states in the shifted BZ
                int n_occupied = 0, n_empty_shifted = 0;
                if (E0 < Ef) n_occupied++;
                if (E1 < Ef) n_occupied++;
                if (E2 < Ef) n_occupied++;
                if (Ep0 > Ef) n_empty_shifted++;
                if (Ep1 > Ef) n_empty_shifted++;
                if (Ep2 > Ef) n_empty_shifted++;

                if (n_occupied > 0 && n_empty_shifted > 0) {
                    double V0 = Ep0 - E0;
                    double V1 = Ep1 - E1;
                    double V2 = Ep2 - E2;

                    bool debug_this = verbose && (total_triangles <= 5);
                    if (debug_this) {
                        std::cout << "Triangle " << total_triangles << ":\n";
                        std::cout << "  k0=" << k0.x << "," << k0.y << " E0=" << E0 << " Ep0=" << Ep0 << "\n";
                        std::cout << "  n_occ=" << n_occupied << " n_empty=" << n_empty_shifted << "\n";
                    }

                    double integral = triangle_integral(V0, V1, V2, triangle_area, debug_this);
                    chi += integral;
                    if (std::abs(integral) > EPSILON) contributing_triangles++;
                }
            }

            // Triangle 2
            {
                total_triangles++;
                Vec3 k0(kx0+dk, ky0, 0);
                Vec3 k1(kx0+dk, ky0+dk, 0);
                Vec3 k2(kx0, ky0+dk, 0);

                double E0 = tb.energy(k0.x, k0.y);
                double E1 = tb.energy(k1.x, k1.y);
                double E2 = tb.energy(k2.x, k2.y);

                double Ep0 = tb.energy(k0.x+q.x, k0.y+q.y);
                double Ep1 = tb.energy(k1.x+q.x, k1.y+q.y);
                double Ep2 = tb.energy(k2.x+q.x, k2.y+q.y);

                int n_occupied = 0, n_empty_shifted = 0;
                if (E0 < Ef) n_occupied++;
                if (E1 < Ef) n_occupied++;
                if (E2 < Ef) n_occupied++;
                if (Ep0 > Ef) n_empty_shifted++;
                if (Ep1 > Ef) n_empty_shifted++;
                if (Ep2 > Ef) n_empty_shifted++;

                if (n_occupied > 0 && n_empty_shifted > 0) {
                    double V0 = Ep0 - E0;
                    double V1 = Ep1 - E1;
                    double V2 = Ep2 - E2;

                    double integral = triangle_integral(V0, V1, V2, triangle_area, false);
                    chi += integral;
                    if (std::abs(integral) > EPSILON) contributing_triangles++;
                }
            }
        }
    }

    chi *= 2.0; // spin factor

    if (verbose) {
        std::cout << "Total triangles: " << total_triangles << "\n";
        std::cout << "Contributing: " << contributing_triangles << "\n";
    }

    return chi;
}

double chi_numerical(const Vec3& q, const TightBinding2D& tb, int Nk, bool verbose=false) {
    double Ef = tb.fermi_energy();
    double chi = 0.0;
    double dk = 2.0*PI / (tb.a * Nk);

    int total_points = 0;
    int contributing_points = 0;

    for (int ix = 0; ix < Nk; ix++) {
        for (int iy = 0; iy < Nk; iy++) {
            total_points++;
            double kx = -PI/tb.a + (ix + 0.5) * dk;
            double ky = -PI/tb.a + (iy + 0.5) * dk;

            double Ek = tb.energy(kx, ky);
            double Ekq = tb.energy(kx+q.x, ky+q.y);

            if (verbose && total_points <= 5) {
                std::cout << "Point " << total_points << ": k=(" << kx << "," << ky
                          << ") E=" << Ek << " Ekq=" << Ekq << "\n";
            }

            // f(Ek)[1-f(Ekq)] at T=0
            if (Ek < Ef && Ekq > Ef) {
                double denom = Ekq - Ek;
                if (std::abs(denom) > EPSILON) {
                    chi += 1.0 / denom;
                    contributing_points++;
                    if (verbose && contributing_points <= 3) {
                        std::cout << "  -> Contributes: 1/" << denom << " = " << 1.0/denom << "\n";
                    }
                }
            }
        }
    }

    chi *= dk * dk * 2.0; // area element * spin

    if (verbose) {
        std::cout << "Total points: " << total_points << "\n";
        std::cout << "Contributing: " << contributing_points << "\n";
    }

    return chi;
}

int main() {
    std::cout << "============================================\n";
    std::cout << "DEBUG: Susceptibility Calculation\n";
    std::cout << "============================================\n\n";

    double mu = -1.0;
    TightBinding2D tb(1.0, 1.0, mu);

    std::cout << "2D Tight-Binding: E(k) = -2t[cos(kx) + cos(ky)] - mu\n";
    std::cout << "  t = " << tb.t << ", a = " << tb.a << ", mu = " << tb.mu << "\n";
    std::cout << "  Ef = " << tb.fermi_energy() << "\n\n";

    int Nk = 40;
    std::cout << "k-mesh: " << Nk << " x " << Nk << "\n\n";

    // Test single q value with verbose output
    std::cout << "========== Testing q = (0.1π, 0.1π) with VERBOSE ==========\n";
    Vec3 q_test(0.1*PI, 0.1*PI, 0);

    std::cout << "\n--- NUMERICAL METHOD ---\n";
    double chi_num_test = chi_numerical(q_test, tb, Nk, true);
    std::cout << "Result: " << chi_num_test << "\n";

    std::cout << "\n--- TETRAHEDRON METHOD ---\n";
    double chi_tet_test = chi_tetrahedron_method(q_test, tb, Nk, true);
    std::cout << "Result: " << chi_tet_test << "\n";

    std::cout << "\n========== Full scan along (π,π) direction ==========\n";
    std::cout << std::setw(12) << "q/(π,π)"
              << std::setw(16) << "χ_tetrahedron"
              << std::setw(16) << "χ_numerical"
              << std::setw(14) << "ratio\n";
    std::cout << std::string(58, '-') << "\n";

    std::ofstream outfile("chi_debug.dat");
    outfile << "# q  chi_tetrahedron  chi_numerical  ratio\n";

    for (int iq = 0; iq <= 20; iq++) {
        double q_mag = (iq * PI) / 20.0;
        Vec3 q(q_mag, q_mag, 0); // along (π,π)

        double chi_tet = chi_tetrahedron_method(q, tb, Nk, false);
        double chi_num = chi_numerical(q, tb, Nk, false);
        double ratio = (std::abs(chi_num) > EPSILON) ? chi_tet/chi_num : 0.0;

        std::cout << std::setw(12) << std::fixed << std::setprecision(4) << q_mag/PI
                  << std::setw(16) << std::setprecision(6) << chi_tet
                  << std::setw(16) << chi_num
                  << std::setw(14) << std::setprecision(4) << ratio << "\n";

        outfile << q_mag/PI << "  " << chi_tet << "  " << chi_num << "  " << ratio << "\n";
    }

    outfile.close();
    std::cout << "\nResults saved to chi_debug.dat\n";

    return 0;
}
