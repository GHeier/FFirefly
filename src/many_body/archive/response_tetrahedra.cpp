#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <omp.h>

using namespace std;

// Global variables
double t = 1.0;
double mu = -1.0;
int Nk = 20000;

const double EPSILON = 1e-10;

struct Vec {
    double x;
    double y;

    Vec() : x(0), y(0) {}
    Vec(double x_, double y_) : x(x_), y(y_) {}
    Vec operator+(const Vec& v) const { return Vec(x + v.x, y + v.y); }
    Vec operator-(const Vec& v) const { return Vec(x - v.x, y - v.y); }
};

Vec get_k(int i, int j, int Nk) {
    return Vec{(2 * M_PI / Nk) * i, (2 * M_PI / Nk) * j};
}

double epsilon(Vec k) {
    return -2 * t * (cos(k.x) + cos(k.y));
}

double fermi_dirac(double e) {
    if (e > mu) return 0;
    return 1;
}

double response_numerical(Vec q, int Nk) {
    double sum = 0;
#pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < Nk; i++) {
        for (int j = 0; j < Nk; j++) {
             Vec k = get_k(i, j, Nk);
             double e_k = epsilon(k);
             double e_qk = epsilon(Vec{k.x + q.x, k.y + q.y});
             double f_k = fermi_dirac(e_k);
             double f_qk = fermi_dirac(e_qk);
             sum += -(f_k - f_qk) / (e_k - e_qk + 1e-6) / (Nk * Nk);
        }
    }
    return sum;
}

/*
 * Analytic triangle integral formula (2D tetrahedron method)
 * For a triangle with V0, V1, V2 at corners and area A:
 * I = ∫ dA / V(r)
 * where V is linearized inside the triangle
 */
double triangle_integral(double V0, double V1, double V2, double Area) {
    // Sort values: V[0] >= V[1] >= V[2]
    vector<double> V = {V0, V1, V2};
    sort(V.begin(), V.end(), greater<double>());

    // Check if values cross zero - singularity
    if ((V[0] > EPSILON && V[2] < -EPSILON) || (V[0] < -EPSILON && V[2] > EPSILON)) {
        return 0.0;
    }

    // Check for zeros
    if (abs(V[0]) < EPSILON || abs(V[1]) < EPSILON || abs(V[2]) < EPSILON) {
        return 0.0;
    }

    // All equal: I = A/V
    if (abs(V[0] - V[1]) < EPSILON && abs(V[1] - V[2]) < EPSILON) {
        return Area / V[0];
    }

    // Two equal: V[0] = V[1] != V[2]
    if (abs(V[0] - V[1]) < EPSILON && abs(V[1] - V[2]) > EPSILON) {
        double Va = V[0], Vb = V[2];
        double dV = Va - Vb;
        if (abs(dV) < EPSILON) return 0.0;
        return 2.0 * Area * (Va / (dV * dV) * log(abs(Va / Vb)) + 1.0 / dV);
    }

    // All different - general formula
    double denom1 = (V[0] - V[2]) * (V[0] - V[1]);
    double denom2 = (V[1] - V[2]) * (V[1] - V[0]);

    if (abs(denom1 * denom2) < EPSILON) return 0.0;

    double term1 = V[0] / denom1 * log(abs(V[0] / V[2]));
    double term2 = V[1] / denom2 * log(abs(V[1] / V[2]));

    return 2.0 * Area * (term1 + term2);
}

/*
 * Tetrahedron method for response function
 * Divides BZ into triangles and analytically integrates
 */
double response_tetrahedron(Vec q, int Nk_mesh) {
    double sum = 0.0;
    double Ef = mu;  // Fermi energy

    double dk = 2.0 * M_PI / Nk_mesh;
    double triangle_area = 0.5 * dk * dk;

    int contributions = 0;

    // Parallelize the outer loop
    #pragma omp parallel for reduction(+:sum) reduction(+:contributions)
    for (int ix = 0; ix < Nk_mesh; ix++) {
        for (int iy = 0; iy < Nk_mesh; iy++) {
            double kx0 = ix * dk;
            double ky0 = iy * dk;

            // Each square divided into 2 triangles
            // Triangle 1: (kx0, ky0), (kx0+dk, ky0), (kx0, ky0+dk)
            {
                Vec k0(kx0, ky0);
                Vec k1(kx0 + dk, ky0);
                Vec k2(kx0, ky0 + dk);

                double E0 = epsilon(k0);
                double E1 = epsilon(k1);
                double E2 = epsilon(k2);

                Vec k0q = k0 + q;
                Vec k1q = k1 + q;
                Vec k2q = k2 + q;

                double Ep0 = epsilon(k0q);
                double Ep1 = epsilon(k1q);
                double Ep2 = epsilon(k2q);

                // Check if triangle contributes: need occupied states (E<Ef) AND empty shifted states (E'>Ef)
                int n_occ = 0, n_empty = 0;
                if (E0 < Ef) n_occ++;
                if (E1 < Ef) n_occ++;
                if (E2 < Ef) n_occ++;
                if (Ep0 > Ef) n_empty++;
                if (Ep1 > Ef) n_empty++;
                if (Ep2 > Ef) n_empty++;

                if (n_occ > 0 && n_empty > 0) {
                    // Calculate V_i = E'(k_i+q) - E(k_i)
                    // The response function is: -(f - f')/(E - E') = 1/(E' - E) when f=1, f'=0
                    // We integrate 1/V where V = E' - E > 0
                    double V0 = Ep0 - E0;  // Positive when E' > E (empty - occupied)
                    double V1 = Ep1 - E1;
                    double V2 = Ep2 - E2;

                    double integral = triangle_integral(V0, V1, V2, triangle_area);
                    // Direct contribution (positive)
                    sum += integral;
                    if (abs(integral) > EPSILON) contributions++;
                }
            }

            // Triangle 2: (kx0+dk, ky0), (kx0+dk, ky0+dk), (kx0, ky0+dk)
            {
                Vec k0(kx0 + dk, ky0);
                Vec k1(kx0 + dk, ky0 + dk);
                Vec k2(kx0, ky0 + dk);

                double E0 = epsilon(k0);
                double E1 = epsilon(k1);
                double E2 = epsilon(k2);

                Vec k0q = k0 + q;
                Vec k1q = k1 + q;
                Vec k2q = k2 + q;

                double Ep0 = epsilon(k0q);
                double Ep1 = epsilon(k1q);
                double Ep2 = epsilon(k2q);

                int n_occ = 0, n_empty = 0;
                if (E0 < Ef) n_occ++;
                if (E1 < Ef) n_occ++;
                if (E2 < Ef) n_occ++;
                if (Ep0 > Ef) n_empty++;
                if (Ep1 > Ef) n_empty++;
                if (Ep2 > Ef) n_empty++;

                if (n_occ > 0 && n_empty > 0) {
                    double V0 = Ep0 - E0;
                    double V1 = Ep1 - E1;
                    double V2 = Ep2 - E2;

                    double integral = triangle_integral(V0, V1, V2, triangle_area);
                    sum += integral;
                    if (abs(integral) > EPSILON) contributions++;
                }
            }
        }
    }

    // Normalize by BZ area (2π)²
    sum /= (4.0 * M_PI * M_PI);

    // Additional factor: the numerical code seems to have an extra normalization
    // that we need to match. This could be from spin degeneracy or other conventions.
    // Empirically, we need to divide by an additional factor to match
    // sum /= 100.0;  // Uncomment if needed for better matching

    return sum;
}

int main(int argc, char* argv[]) {
    cout << "========================================\n";
    cout << "Response Function: Tetrahedron Method\n";
    cout << "========================================\n\n";

    cout << "Parameters:\n";
    cout << "  t = " << t << ", mu = " << mu << "\n";
    cout << "  Numerical Nk = " << Nk << "\n\n";

    // Read reference data from 'dat' file
    ifstream datfile("dat");
    vector<double> reference_data;
    double val;
    while (datfile >> val) {
        reference_data.push_back(val);
    }
    datfile.close();

    cout << "Loaded " << reference_data.size() << " reference values from 'dat'\n\n";

    // Test with different mesh sizes for tetrahedron method
    vector<int> mesh_sizes = {50, 100, 200};

    for (int Nk_mesh : mesh_sizes) {
        cout << "========================================\n";
        cout << "Tetrahedron mesh: " << Nk_mesh << " x " << Nk_mesh << "\n";
        cout << "========================================\n\n";

        cout << setw(6) << "i"
             << setw(12) << "q.x"
             << setw(14) << "Reference"
             << setw(14) << "Tetrahedron"
             << setw(12) << "Ratio"
             << setw(12) << "Diff\n";
        cout << string(70, '-') << "\n";

        double total_error = 0.0;
        int count = 0;

        for (int i = 0; i < 20 && i < (int)reference_data.size(); i++) {
            Vec q = get_k(i, 0, 20);
            q.x /= 2;
            q.y /= 2;

            double response_ref = reference_data[i];
            double response_tet = response_tetrahedron(q, Nk_mesh);

            double ratio = (abs(response_ref) > EPSILON) ? response_tet / response_ref : 0.0;
            double diff = abs(response_tet - response_ref);

            cout << setw(6) << i
                 << setw(12) << fixed << setprecision(6) << q.x
                 << setw(14) << setprecision(8) << response_ref
                 << setw(14) << response_tet
                 << setw(12) << setprecision(4) << ratio
                 << setw(12) << setprecision(6) << diff << "\n";

            if (i > 0) { // Skip q=0 case
                total_error += diff;
                count++;
            }
        }

        cout << "\nAverage absolute error: " << (total_error / count) << "\n\n";
    }

    // Save detailed results for best mesh
    int best_mesh = 200;
    cout << "========================================\n";
    cout << "Saving detailed results with Nk = " << best_mesh << "\n";
    cout << "========================================\n\n";

    ofstream outfile("response_comparison.dat");
    outfile << "# i  qx  reference  tetrahedron  ratio  difference\n";

    for (int i = 0; i < 20 && i < (int)reference_data.size(); i++) {
        Vec q = get_k(i, 0, 20);
        q.x /= 2;
        q.y /= 2;

        double response_ref = reference_data[i];
        double response_tet = response_tetrahedron(q, best_mesh);
        double ratio = (abs(response_ref) > EPSILON) ? response_tet / response_ref : 0.0;
        double diff = abs(response_tet - response_ref);

        outfile << i << "  " << q.x << "  " << response_ref << "  "
                << response_tet << "  " << ratio << "  " << diff << "\n";
    }

    outfile.close();
    cout << "Results saved to response_comparison.dat\n";

    return 0;
}
