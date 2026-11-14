#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

#include "src/algorithms/fft.hpp"
#include "src/config/load/cpp_config.hpp"
#include "src/objects/CMField/fields.hpp"
#include "src/objects/CMField/hamiltonian.hpp"
#include "src/objects/vec.hpp"
#include "self_energy.hpp"

using namespace std;

static void print_progress(int current, int total, std::chrono::steady_clock::time_point start) {
    double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    double rate = current / std::max(elapsed, 1e-9);
    double remaining = (total - current) / std::max(rate, 1e-9);
    if (remaining > 60.0)
        printf("  Progress: %d/%d | Rate: %.1f/s | ETA: %.1f min        \r", current, total, rate, remaining / 60.0);
    else
        printf("  Progress: %d/%d | Rate: %.1f/s | ETA: %.1f s           \r", current, total, rate, remaining);
    fflush(stdout);
}

/**
 * T=0 Kramers-Kronig transform: recovers Re Sigma(q, w_i) from Im Sigma(q, w_j)
 * tabulated on the same (possibly non-uniform) frequency grid, via
 *   Re Sigma(w) = (1/pi) P integral dw' Im Sigma(w') / (w' - w)
 * The (w_j - w_i) -> 0 singularity is regularized by subtracting Im Sigma(w_i)
 * from the numerator before summing; this is exact for the principal value as
 * long as Im Sigma has decayed to ~0 near the edges of the grid (i.e. the
 * cutoff_energy/w_pts used to build the vertex file is large/dense enough).
 */
static vector<float> kk_transform(const vector<float>& imS, const vector<float>& wgrid) {
    int n = (int)wgrid.size();
    vector<float> reS(n, 0.0f);
    if (n < 2) return reS;

    vector<float> dw(n);
    dw[0] = wgrid[1] - wgrid[0];
    dw[n - 1] = wgrid[n - 1] - wgrid[n - 2];
    for (int j = 1; j < n - 1; j++) dw[j] = 0.5f * (wgrid[j + 1] - wgrid[j - 1]);

    for (int i = 0; i < n; i++) {
        float acc = 0.0f;
        for (int j = 0; j < n; j++) {
            if (j == i) continue;
            acc += (imS[j] - imS[i]) * dw[j] / (wgrid[j] - wgrid[i]);
        }
        reS[i] = acc / (float)M_PI;
    }
    return reS;
}

/**
 * Trapezoidal integration weights for a (possibly non-uniform) 1D grid.
 */
static vector<float> trapezoidal_weights(const vector<float>& grid) {
    int n = (int)grid.size();
    vector<float> dw(n, 0.0f);
    if (n < 2) return dw;
    dw[0] = grid[1] - grid[0];
    dw[n - 1] = grid[n - 1] - grid[n - 2];
    for (int j = 1; j < n - 1; j++) dw[j] = 0.5f * (grid[j + 1] - grid[j - 1]);
    return dw;
}

/**
 * Fermi and Bose occupation functions. T<=0 reproduces the exact T=0 step-
 * function limit (Theta(mu-eps), -Theta(-v)) bit-for-bit, so callers can use
 * these unconditionally without a separate T=0 code path. Callers must
 * handle bose()'s v=0 removable singularity themselves (see the L'Hopital
 * handling in call_flex_self_energy()/call_flex_self_energy_fft()) - bose()
 * does not special-case v~=0 on its own.
 */
static inline float fermi(float eps, float mu, float T) {
    if (T <= 0.0f) return (eps < mu) ? 1.0f : 0.0f;
    float x = (eps - mu) / T;
    if (x > 40.0f) return 0.0f;
    if (x < -40.0f) return 1.0f;
    return 1.0f / (std::exp(x) + 1.0f);
}

static inline float bose(float v, float T) {
    if (T <= 0.0f) return (v < 0.0f) ? -1.0f : 0.0f;
    float x = v / T;
    if (x > 40.0f) return 0.0f;
    if (x < -40.0f) return -1.0f;
    return 1.0f / (std::exp(x) - 1.0f);
}

/**
 * FFT-based Kramers-Kronig transform: recovers Re Sigma(w_i) from Im Sigma(w_j)
 * via the standard "analytic signal" sign-filter method,
 *   Re Sigma = IFFT[ i*sgn(xi) * FFT(Im Sigma) ],
 * rather than literally transforming the singular kernel h(w)=1/(pi*w)
 * (numerically fragile at its own singularity). Requires wgrid to be
 * uniformly spaced for the FFT convolution to represent the correct
 * continuous integral - falls back to the exact kk_transform() otherwise.
 * Zero-pads to the next power of two >= 2n to avoid circular wraparound,
 * since the KK relation is not a periodic convolution.
 *
 * O(n log n) vs kk_transform()'s O(n^2). The exact sign/normalization
 * convention here was validated empirically against kk_transform() on real
 * data (see from_vertex README) rather than trusted from derivation alone -
 * kk_transform() is kept as the correctness reference.
 *
 * Thread safety: fft_1d_f() creates and destroys an FFTW plan on every call
 * (see fft.cpp), and FFTW's planner is NOT thread-safe by default (this repo
 * doesn't link fftw3f_threads/fftw3f_omp, which is what would let
 * fftwf_make_planner_thread_safe() be used instead) - calling it
 * concurrently from multiple OpenMP threads (this function is called once
 * per q-point inside call_flex_self_energy()'s parallelized q-loop)
 * corrupts FFTW's internal planner state (observed directly: "fftw:
 * planner.c:261: assertion failed" followed by a segfault). The two FFT
 * calls below are therefore serialized via a named critical section - the
 * O(n log n) work itself is otherwise unaffected, so this only serializes
 * plan creation/execution, not the (already fast) surrounding computation.
 */
static vector<float> kk_transform_fft(const vector<float>& imS, const vector<float>& wgrid) {
    int n = (int)wgrid.size();
    if (n < 2) return vector<float>(n, 0.0f);

    float dw0 = wgrid[1] - wgrid[0];
    for (int i = 1; i < n - 1; i++) {
        if (std::abs((wgrid[i + 1] - wgrid[i]) - dw0) > 1e-4f * std::abs(dw0))
            return kk_transform(imS, wgrid);
    }

    int m = 1;
    while (m < 2 * n) m *= 2;

    vector<complex<float>> data(m, complex<float>(0.0f, 0.0f));
    for (int i = 0; i < n; i++) data[i] = complex<float>(imS[i], 0.0f);

    vector<complex<float>> F;
    #pragma omp critical(fftw_planner)
    F = fft_1d_f(data, false, false);

    for (int k = 0; k < m; k++) {
        if (k == 0 || (m % 2 == 0 && k == m / 2)) {
            F[k] = complex<float>(0.0f, 0.0f);
        } else if (k < m / 2) {
            F[k] *= complex<float>(0.0f, 1.0f);
        } else {
            F[k] *= complex<float>(0.0f, -1.0f);
        }
    }

    vector<complex<float>> result;
    #pragma omp critical(fftw_planner)
    result = fft_1d_f(F, true, false);

    vector<float> reS(n);
    for (int i = 0; i < n; i++) reS[i] = result[i].real();
    return reS;
}

/**
 * Shared tail: given Sigma(q,w) already computed (by either the direct sum
 * or the FFT convolution), print diagnostics, save it, build the retarded
 * Green's function via the Dyson equation, and save that too.
 */
static void finish_self_energy(vector<vector<cfloat>>& sigma, const vector<Vec>& qpts,
                                const vector<float>& wpts, int chidim, bool q_as_mesh,
                                const vector<int>& mesh_dims,
                                const vector<vector<float>>& stored_points, Hamiltonian& H) {
    int Nq = (int)qpts.size();
    int nw = (int)wpts.size();

    int n_nonfinite = 0;
    for (int iq = 0; iq < Nq; iq++)
        for (int iw = 0; iw < nw; iw++)
            if (!std::isfinite(sigma[iq][iw].real()) || !std::isfinite(sigma[iq][iw].imag()))
                n_nonfinite++;
    if (n_nonfinite > 0)
        printf("WARNING: %d/%d Sigma entries are NaN/Inf - check the v~=0 L'Hopital handling"
               " (Temperature>0) or the input vertex file.\n", n_nonfinite, Nq * nw);

    float max_im = 0.0f, max_re = 0.0f;
    for (int iq = 0; iq < Nq; iq++) {
        for (int iw = 0; iw < nw; iw++) {
            if (!std::isfinite(sigma[iq][iw].imag()) || !std::isfinite(sigma[iq][iw].real())) continue;
            max_im = std::max(max_im, std::abs(sigma[iq][iw].imag()));
            max_re = std::max(max_re, std::abs(sigma[iq][iw].real()));
        }
    }
    printf("Max |Im Sigma|: %f\n", max_im);
    printf("Max |Re Sigma|: %f\n", max_re);

    // Flatten to (w outer, q inner) order, matching the w-k convention
    // already used for chi/vertex mesh files.
    vector<cfloat> sigma_vals;
    sigma_vals.reserve((size_t)Nq * nw);
    for (int iw = 0; iw < nw; iw++)
        for (int iq = 0; iq < Nq; iq++)
            sigma_vals.push_back(sigma[iq][iw]);

    printf("Saving Self-Energy\n");
    string file = outdir + prefix + "_self_energy." + filetype;
    if (filetype == "hdf5" || filetype == "h5") {
        if (q_as_mesh) {
            save_data(file, sigma_vals, vector<int>{}, mesh_dims, brillouin_zone, wpts);
        } else {
            save_data(file, sigma_vals, vector<int>{}, vector<int>{}, vector<vector<float>>{{}}, wpts, stored_points);
        }
    }
    cout << "Saved to " << file << endl;

    // Retarded Green's function via the Dyson equation, built directly from
    // the complex self-energy just computed above:
    //   G(q,w) = 1 / (w - e(q) - Sigma(q,w))
    // No extra +i*eta broadening is added - causality is already carried by
    // Im Sigma(q,w) <= 0. At a mesh point where Im Sigma(q,w) happens to be
    // exactly 0 and Re[w-e(q)-Sigma] also vanishes, G diverges - a genuine
    // undamped-pole artifact, not a numerical error.
    printf("Computing G(q,w) = 1/(w - e(q) - Sigma(q,w))\n");
    vector<vector<cfloat>> Gf(Nq, vector<cfloat>(nw));
    #pragma omp parallel for schedule(dynamic)
    for (int iq = 0; iq < Nq; iq++) {
        float e_q = H.get_bands(qpts[iq])[0];
        for (int iw = 0; iw < nw; iw++) {
            cfloat denom = cfloat(wpts[iw] - e_q, 0.0f) - sigma[iq][iw];
            Gf[iq][iw] = cfloat(1.0f, 0.0f) / denom;
        }
    }

    float max_G = 0.0f;
    for (int iq = 0; iq < Nq; iq++) {
        for (int iw = 0; iw < nw; iw++) {
            if (std::isfinite(Gf[iq][iw].real()) && std::isfinite(Gf[iq][iw].imag()))
                max_G = std::max(max_G, std::abs(Gf[iq][iw]));
        }
    }
    printf("Max |G| (finite entries only): %f\n", max_G);

    vector<cfloat> G_vals;
    G_vals.reserve((size_t)Nq * nw);
    for (int iw = 0; iw < nw; iw++)
        for (int iq = 0; iq < Nq; iq++)
            G_vals.push_back(Gf[iq][iw]);

    printf("Saving Green's function\n");
    string Gfile = outdir + prefix + "_G." + filetype;
    if (filetype == "hdf5" || filetype == "h5") {
        if (q_as_mesh) {
            save_data(Gfile, G_vals, vector<int>{}, mesh_dims, brillouin_zone, wpts);
        } else {
            save_data(Gfile, G_vals, vector<int>{}, vector<int>{}, vector<vector<float>>{{}}, wpts, stored_points);
        }
    }
    cout << "Saved to " << Gfile << endl;

    // Spectral function A(q,w) = -2 Im G(q,w). Real and >= 0 for a physical
    // retarded G (equivalent to Im Sigma(q,w) <= 0); a negative value here
    // signals either a causality violation in Sigma or an interpolation
    // artifact when this file is later evaluated off the native mesh (linear
    // interpolation of a near-singular complex function doesn't preserve the
    // causal 1/(a-ib) structure near a pole - see the from_vertex README).
    printf("Computing A(q,w) = -2 Im G(q,w)\n");
    vector<float> A_vals;
    A_vals.reserve((size_t)Nq * nw);
    float max_A = 0.0f, min_A = 0.0f;
    for (int iq = 0; iq < Nq; iq++) {
        for (int iw = 0; iw < nw; iw++) {
            float a = -2.0f * Gf[iq][iw].imag();
            if (std::isfinite(a)) {
                max_A = std::max(max_A, a);
                min_A = std::min(min_A, a);
            }
        }
    }
    for (int iw = 0; iw < nw; iw++)
        for (int iq = 0; iq < Nq; iq++)
            A_vals.push_back(-2.0f * Gf[iq][iw].imag());

    printf("Max A: %f, Min A: %f%s\n", max_A, min_A,
           (min_A < -1e-4f) ? "  (negative - see causality note above)" : "");

    printf("Saving spectral function\n");
    string Afile = outdir + prefix + "_A_w." + filetype;
    if (filetype == "hdf5" || filetype == "h5") {
        if (q_as_mesh) {
            save_data(Afile, A_vals, vector<int>{}, mesh_dims, brillouin_zone, wpts);
        } else {
            save_data(Afile, A_vals, vector<int>{}, vector<int>{}, vector<vector<float>>{{}}, wpts, stored_points);
        }
    }
    cout << "Saved to " << Afile << endl;
}

/**
 * Real-axis FLEX/spin-fluctuation self-energy from the vertex V(k,w), at
 * temperature Temperature (global config extern; Temperature<=0 is the exact
 * T=0 limit):
 *
 *   Sigma(q,w) = (1/N) sum_k integral (dv/pi) ImV(k,v) [n_B(v)+f(e(q-k))] / (w-v-e(q-k)+i*eta)
 *
 * The v-integral collapses analytically onto v = w - e(q-k) via the delta-
 * function identity Im[1/(x+i*eta)] -> -pi*delta(x) (no broadening/
 * tetrahedron resolution needed - V(k,v) is evaluated directly via its own
 * frequency interpolation, at ANY temperature, not just T=0), giving:
 *
 *   Im Sigma(q,w)  = -(1/N) sum_k [n_B(w-e(q-k)) + f(e(q-k))] * ImV(k, w-e(q-k))
 *   Re Sigma(q,w)  = Kramers-Kronig transform of Im Sigma(q,w) over w
 *
 * At Temperature<=0, n_B(v) -> -Theta(-v) and f(e) -> Theta(mu-e) exactly
 * (see fermi()/bose()), reproducing the original T=0 step-function formula
 * bit-for-bit. Re Sigma is obtained via kk_transform_fft() (an O(n log n)
 * FFT-based Hilbert transform, falling back to the exact O(n^2)
 * kk_transform() for a non-uniform w-grid) - this is temperature-agnostic,
 * only consuming the already-computed Im Sigma.
 *
 * The retarded Green's function is then built from Sigma via the Dyson
 * equation, G(q,w) = 1 / (w - e(q) - Sigma(q,w)); the spectral function
 * A(q,w) = -2 Im G(q,w) is computed directly from G. Both are saved
 * alongside Sigma (see finish_self_energy()).
 *
 * q_mesh sets the output grid Sigma/G/A are evaluated on; k_mesh sets the
 * internal BZ-sum grid above (where V is evaluated - "V(k)"). The two are
 * independent - changing q_mesh changes only the resolution of the saved
 * Sigma(q,w)/G(q,w), while changing k_mesh changes the actual physics
 * (BZ-integration accuracy of the sum).
 */
void call_flex_self_energy() {
    string filename = outdir + prefix + "_vertex_singlet." + filetype;
    printf("Reading vertex from %s\n", filename.c_str());

    Field_C V(filename);
    Hamiltonian H;
    float mu = fermi_energy;
    int chidim = V.get_data()->dimension;

    vector<float> wpts = V.get_data()->w_points;
    if (wpts.size() == 0) {
        wpts.push_back(0.0f);
    }
    int nw = (int)wpts.size();

    printf("Self-energy (direct sum): T = %f\n", Temperature);
    printf("  mu = %f\n", mu);

    // Internal k-mesh used for the BZ sum - always a full mesh over the
    // configured k_mesh/brillouin_zone, independent of how the output
    // q-points below are chosen. V(k,w) is evaluated at these points.
    float nkx = (float)k_mesh[0], nky = (float)k_mesh[1], nkz = (float)k_mesh[2];
    if (chidim == 2) nkz = 1;
    vector<Vec> kpts;
    kpts.reserve((size_t)(nkx * nky * nkz));
    for (int i = 0; i < (int)nkx; i++) {
        for (int j = 0; j < (int)nky; j++) {
            for (int k = 0; k < (int)nkz; k++) {
                Vec kp = brillouin_zone * Vec(i / nkx - 0.5f, j / nky - 0.5f, k / nkz - 0.5f);
                kp.dimension = chidim;
                kpts.push_back(kp);
            }
        }
    }
    int Nk = (int)kpts.size();

    // Output q-points: a mesh built from q_mesh, independent of the k_mesh
    // used for the internal BZ sum above - Sigma/G resolution and BZ-sum
    // accuracy are separate knobs (q_mesh controls the former, k_mesh the
    // latter). V(k,w) is evaluated via its own interpolation, so the output
    // q-mesh doesn't need to match whatever mesh V/chi were computed on.
    // Falls back to V's own stored points only when V itself isn't on a
    // regular mesh (as_mesh=false), since a fresh q-mesh wouldn't have any
    // defined correspondence to scattered data.
    bool q_as_mesh = V.get_data()->as_mesh;
    vector<Vec> qpts;
    vector<int> mesh_dims;
    if (q_as_mesh) {
        float nqx = (float)q_mesh[0], nqy = (float)q_mesh[1], nqz = (float)q_mesh[2];
        if (chidim == 2) nqz = 1;
        qpts.reserve((size_t)(nqx * nqy * nqz));
        for (int i = 0; i < (int)nqx; i++) {
            for (int j = 0; j < (int)nqy; j++) {
                for (int k = 0; k < (int)nqz; k++) {
                    Vec qp = brillouin_zone * Vec(i / nqx - 0.5f, j / nqy - 0.5f, k / nqz - 0.5f);
                    qp.dimension = chidim;
                    qpts.push_back(qp);
                }
            }
        }
        mesh_dims.assign(q_mesh.begin(), q_mesh.begin() + chidim);
    } else {
        printf("Using stored point data for output q-points\n");
        for (const auto& point : V.get_data()->points) {
            Vec q(point.data(), chidim);
            q.dimension = chidim;
            qpts.push_back(q);
        }
    }
    int Nq = (int)qpts.size();
    printf("Computing Sigma(q,w) on %d q-points x %d frequencies (k-sum over %d points)\n", Nq, nw, Nk);

    // band index 0 only (single/diagonal-band approximation, matching the
    // rest of the chi0/vertex pipeline - extend with a band loop if nbnd>1
    // matters for a given material).
    vector<vector<cfloat>> sigma(Nq, vector<cfloat>(nw));

    // V(k,w) is only tabulated over [w_min, w_max] (the cutoff_energy window
    // used to build the vertex/chi files); the field evaluator does not
    // extrapolate past that. Physically Im V(k,nu) is negligible once |nu|
    // exceeds the particle-hole bandwidth captured by that window, so any
    // nu = w - e(q-k) landing outside it is simply dropped from the sum
    // rather than extrapolated.
    float w_min = wpts.front(), w_max = wpts.back();
    printf("V(k,w) tabulated on [%f, %f]; terms with |w - e(q-k)| outside this window are dropped\n", w_min, w_max);

    // Precompute d(ImV)/dv at v=0 for each k, needed only at T>0 to handle
    // bose(v,T)'s removable singularity at v=0 via L'Hopital's rule: ImV(k,0)=0
    // for all k by causality (Im V is odd in frequency), so bose(v,T)*ImV(k,v)
    // is a 0/0 form as v->0 whose correct value is -T * d(ImV(k,v))/dv|_0.
    // Depends only on k, so computed once here rather than per (q,w).
    vector<float> dImV_dv0;
    float v_eps = 0.0f;
    if (Temperature > 0.0f) {
        v_eps = 1e-3f * Temperature;
        float delta = 0.01f * (w_max - w_min) / std::max(1, nw - 1);
        dImV_dv0.resize(Nk);
        #pragma omp parallel for
        for (int ik = 0; ik < Nk; ik++) {
            dImV_dv0[ik] = (V(kpts[ik], delta).imag() - V(kpts[ik], -delta).imag()) / (2.0f * delta);
        }
    }

    // Parallelism lives on the outer q-loop (one thread team for the whole
    // loop, dynamic-scheduled) - measured ~20-30% faster than parallelizing
    // the inner k-sum, which pays repeated fork/join/barrier overhead once
    // per (q,w) pair. Progress is still reported via a shared atomic counter
    // (same pattern as response/tetrahedra/response_utils.jl's
    // print_progress), printed by whichever thread crosses a reporting
    // threshold, guarded by a critical section.
    std::atomic<int> progress_counter{0};
    int report_every = std::max(1, Nq / 200);
    auto start_time = std::chrono::steady_clock::now();

    #pragma omp parallel for schedule(dynamic)
    for (int iq = 0; iq < Nq; iq++) {
        Vec q = qpts[iq];
        vector<float> imS(nw, 0.0f);

        for (int iw = 0; iw < nw; iw++) {
            float w = wpts[iw];
            float sum = 0.0f;
            for (int ik = 0; ik < Nk; ik++) {
                Vec qk = q - kpts[ik];
                qk.dimension = chidim;
                float e_qk = H.get_bands(qk)[0];

                float nu = w - e_qk;
                if (nu < w_min || nu > w_max) continue;

                if (Temperature > 0.0f && std::abs(nu) < v_eps) {
                    // bose(nu,Temperature)*ImV(k,nu) -> L'Hopital limit as nu->0
                    // (see dImV_dv0 comment above); fermi() has no singularity here.
                    sum += fermi(e_qk, mu, Temperature) * V(kpts[ik], nu).imag()
                         - Temperature * dImV_dv0[ik];
                    continue;
                }

                float occ = fermi(e_qk, mu, Temperature) + bose(nu, Temperature);
                if (std::abs(occ) < 1e-6f) continue;  // exponential-tail cutoff (perf only)

                cfloat Vknu = V(kpts[ik], nu);
                sum += occ * Vknu.imag();
            }
            imS[iw] = -sum / Nk;
        }

        vector<float> reS = kk_transform_fft(imS, wpts);
        for (int iw = 0; iw < nw; iw++)
            sigma[iq][iw] = cfloat(reS[iw], imS[iw]);

        int done = ++progress_counter;
        if (done % report_every == 0 || done == Nq) {
            #pragma omp critical
            print_progress(done, Nq, start_time);
        }
    }
    printf("\n");

    finish_self_energy(sigma, qpts, wpts, chidim, q_as_mesh, mesh_dims, V.get_data()->points, H);
}

/**
 * FFT-accelerated variant of call_flex_self_energy(), using the *un-collapsed*
 * general formula
 *
 *   Sigma(q,w_i) = (1/N) sum_k sum_j (dv_j/pi) ImV(k,v_j) [n_B(v_j)+f(e(q-k))] / (w_i-v_j-e(q-k)+i*eta)
 *
 * For each fixed pair (v_j, w_i), the occupation factor n_B(v_j)+f(e(q-k))
 * depends on r=q-k only through f(e(r)) - n_B(v_j) is a constant for fixed
 * v_j, at ANY temperature (not just T=0, where n_B(v_j)->-Theta(-v_j) and
 * f(e)->Theta(mu-e) exactly - see fermi()/bose()) - so the k-sum
 *
 *   sum_k ImV(k,v_j) * D_{j,i}(q-k),   D_{j,i}(r) = [f(e(r))+n_B(v_j)] / (w_i-v_j-e(r)+i*eta)
 *
 * is a genuine cyclic convolution over the BZ mesh, computable via
 * IFFT(FFT(ImV(.,v_j)) * FFT(D_{j,i})) in O(Nk log Nk) instead of O(Nk) per
 * q-point (O(Nk^2) total). Summed over the nw x nw grid of (v_j,w_i) pairs,
 * total cost is O(Nw^2 * Nk log Nk) instead of the direct sum's
 * O(Nq*Nk*Nw) - a large win once Nk is big, since Nw is typically small
 * (~10-20 frequency points) while Nk can be 10^4-10^6.
 *
 * Trade-off: unlike the direct method's exact analytic delta-function
 * resolution, this needs a real broadening eta (a genuine Lorentzian, not
 * just a numerical nudge) to keep 1/(w_i-v_j-e(r)+i*eta) well-defined. eta
 * must be resolved by the k-mesh's own energy spacing (bandwidth/Nk^(1/dim))
 * or the discrete sum misses the Lorentzian almost everywhere except rare
 * near-coincidences, giving noise instead of a controlled approximation -
 * so eta is set adaptively from the mesh rather than a fixed constant.
 * Results converge to the direct method's as the mesh is refined (Nk -> inf,
 * eta -> 0 together) and as V's own w_points grid gets denser, but are not
 * bit-identical to it at finite resolution.
 *
 * Requires k_mesh == q_mesh: the q-k index arithmetic that FFT convolution
 * performs implicitly only represents true BZ momentum subtraction when k
 * and q live on the exact same mesh (same dims, same brillouin_zone, same
 * point ordering). q_mesh is therefore ignored in this path - k_mesh sets
 * both the internal sum grid and the output grid. Also requires V to be a
 * regular mesh field (as_mesh=true); FFT convolution has no meaning against
 * scattered stored points.
 */
void call_flex_self_energy_fft() {
    string filename = outdir + prefix + "_vertex_singlet." + filetype;
    printf("Reading vertex from %s\n", filename.c_str());

    Field_C V(filename);
    Hamiltonian H;
    float mu = fermi_energy;
    int chidim = V.get_data()->dimension;

    if (!V.get_data()->as_mesh) {
        printf("ERROR: call_flex_self_energy_fft() requires V to be on a regular mesh (as_mesh=true).\n");
        return;
    }

    vector<float> wpts = V.get_data()->w_points;
    if (wpts.size() == 0) {
        wpts.push_back(0.0f);
    }
    int nw = (int)wpts.size();
    vector<float> dnu = trapezoidal_weights(wpts);

    printf("Self-energy (FFT convolution): T = %f\n", Temperature);
    printf("  mu = %f\n", mu);
    printf("  Note: q_mesh is ignored in FFT mode - k_mesh sets both the internal sum grid\n");
    printf("  and the output grid, since FFT convolution requires them to coincide.\n");

    // Single shared mesh for both the k-sum and the output q-grid.
    int nx = k_mesh[0], ny = k_mesh[1], nz = (chidim == 2) ? 1 : k_mesh[2];
    vector<int> dims = (chidim == 2) ? vector<int>{nx, ny} : vector<int>{nx, ny, nz};
    vector<Vec> kpts;
    kpts.reserve((size_t)nx * ny * nz);
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                Vec kp = brillouin_zone * Vec((float)i / nx - 0.5f, (float)j / ny - 0.5f, (float)k / nz - 0.5f);
                kp.dimension = chidim;
                kpts.push_back(kp);
            }
        }
    }
    int Nk = (int)kpts.size();
    const vector<Vec>& qpts = kpts;  // same mesh, same ordering
    int Nq = Nk;
    vector<int> mesh_dims(k_mesh.begin(), k_mesh.begin() + chidim);
    printf("Computing Sigma(q,w) on %d q-points x %d frequencies (FFT k-sum over %d points)\n", Nq, nw, Nk);

    // Band energy at every mesh point (shared by both k and q roles, same mesh).
    vector<float> eps(Nk);
    #pragma omp parallel for
    for (int i = 0; i < Nk; i++) eps[i] = H.get_bands(kpts[i])[0];

    // Adaptive Lorentzian broadening: eta must be resolved by the k-mesh's
    // own energy spacing, or the discrete sum over mesh points effectively
    // misses the Lorentzian everywhere except the rare near-coincidences,
    // giving noise dominated by mesh-alignment luck rather than a controlled
    // approximation to the delta function. Typical spacing between adjacent
    // mesh points in energy scales as bandwidth / Nk^(1/dim), so eta is set
    // to a fraction of that - shrinks automatically as the mesh is refined
    // (Nk -> infinity, eta -> 0, converging to the direct method's exact
    // result), never smaller than a small floor for degenerate/flat bands.
    float eps_min = *std::min_element(eps.begin(), eps.end());
    float eps_max = *std::max_element(eps.begin(), eps.end());
    float bandwidth = eps_max - eps_min;
    float eta = std::max(1e-3f, bandwidth / (2.0f * std::pow((float)Nk, 1.0f / chidim)));
    printf("  eta = %f (adaptive: bandwidth=%f, Nk=%d)\n", eta, bandwidth, Nk);

    // Im V(k, v_j) on the mesh, for every tabulated v_j - read once via V's
    // own interpolation (exact at these grid points), reused across all
    // (v_j, w_i) pairs below.
    vector<vector<float>> ImV(nw, vector<float>(Nk));
    for (int j = 0; j < nw; j++)
        for (int i = 0; i < Nk; i++)
            ImV[j][i] = V(kpts[i], wpts[j]).imag();

    // FFT of ImV(.,v_j) for each v_j, computed once and reused for every w_i.
    vector<vector<complex<float>>> ImV_fft(nw, vector<complex<float>>(Nk));
    for (int j = 0; j < nw; j++) {
        vector<complex<float>> data(Nk);
        for (int i = 0; i < Nk; i++) data[i] = complex<float>(ImV[j][i], 0.0f);
        ImV_fft[j] = fft_nd_f(data, dims, false);
    }

    // v~=0 special case: bose(v_j,T)'s removable singularity at v_j=0 (see
    // dImV_dv0 comment in call_flex_self_energy()) needs L'Hopital handling
    // when V's own frequency grid includes (or nearly includes) v=0 - common,
    // since w-point grids are often built symmetric around 0. Precompute
    // -T * d(ImV(.,v))/dv|_0 from the already-tabulated ImV table (no new
    // V() evaluations needed) and its FFT once, reused for every iw below.
    int j0 = -1;
    vector<complex<float>> dImV_dv0_fft;
    if (Temperature > 0.0f) {
        float v_eps = 1e-3f * Temperature;
        int best = 0;
        for (int j = 1; j < nw; j++)
            if (std::abs(wpts[j]) < std::abs(wpts[best])) best = j;
        if (std::abs(wpts[best]) < v_eps && best > 0 && best < nw - 1) {
            j0 = best;
            vector<complex<float>> dImV_dv0(Nk);
            float dw = wpts[j0 + 1] - wpts[j0 - 1];
            for (int r = 0; r < Nk; r++)
                dImV_dv0[r] = complex<float>(-Temperature * (ImV[j0 + 1][r] - ImV[j0 - 1][r]) / dw, 0.0f);
            dImV_dv0_fft = fft_nd_f(dImV_dv0, dims, false);
        }
    }

    vector<vector<cfloat>> sigma(Nq, vector<cfloat>(nw));

    auto start_time = std::chrono::steady_clock::now();
    for (int iw = 0; iw < nw; iw++) {
        float w = wpts[iw];
        vector<complex<float>> accum(Nk, complex<float>(0.0f, 0.0f));

        for (int j = 0; j < nw; j++) {
            float nu = wpts[j];
            // At j==j0, the divergent bose(v_j) piece is dropped here and
            // added back separately below via the L'Hopital replacement.
            float nB_j = (j == j0) ? 0.0f : bose(nu, Temperature);

            // D_{j,i}(r) = [f(e(r)) + n_B(v_j)] / (w - v_j - e(r) + i*eta)
            vector<complex<float>> D(Nk);
            for (int r = 0; r < Nk; r++) {
                float occ = fermi(eps[r], mu, Temperature) + nB_j;
                complex<float> denom(w - nu - eps[r], eta);
                D[r] = occ / denom;
            }
            vector<complex<float>> D_fft = fft_nd_f(D, dims, false);

            // Pointwise product in transformed space, then inverse FFT to
            // recover the cyclic convolution sum_k ImV(k,v_j)*D_{j,i}(q-k).
            vector<complex<float>> product(Nk);
            for (int r = 0; r < Nk; r++) product[r] = ImV_fft[j][r] * D_fft[r];
            vector<complex<float>> conv = fft_nd_f(product, dims, true);

            float weight = dnu[j] / (float)M_PI;
            for (int q = 0; q < Nk; q++) accum[q] += weight * conv[q];

            if (j == j0) {
                // L'Hopital replacement for the dropped bose(v_j)*ImV(k,v_j)
                // term: convolve -T*d(ImV)/dv|_0 against 1/denom only (no
                // occupation factor - see dImV_dv0_fft comment above).
                vector<complex<float>> Dinv(Nk);
                for (int r = 0; r < Nk; r++)
                    Dinv[r] = complex<float>(1.0f, 0.0f) / complex<float>(w - nu - eps[r], eta);
                vector<complex<float>> Dinv_fft = fft_nd_f(Dinv, dims, false);
                vector<complex<float>> product2(Nk);
                for (int r = 0; r < Nk; r++) product2[r] = dImV_dv0_fft[r] * Dinv_fft[r];
                vector<complex<float>> conv2 = fft_nd_f(product2, dims, true);
                for (int q = 0; q < Nk; q++) accum[q] += weight * conv2[q];
            }
        }

        for (int q = 0; q < Nq; q++)
            sigma[q][iw] = cfloat(accum[q].real() / Nk, accum[q].imag() / Nk);

        print_progress(iw + 1, nw, start_time);
    }
    printf("\n");

    finish_self_energy(sigma, qpts, wpts, chidim, true, mesh_dims, V.get_data()->points, H);
}
