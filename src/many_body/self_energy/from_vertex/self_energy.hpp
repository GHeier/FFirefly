#pragma once

// Direct real-axis sum at temperature Temperature (global config extern;
// Temperature<=0 is the exact T=0 limit), exact analytic delta-function
// resolution of the frequency integral (no broadening): O(Nq*Nk*Nw) for the
// main sum, plus O(Nq*Nw*log(Nw)) for the FFT-based Kramers-Kronig transform
// that recovers Re Sigma from Im Sigma.
void call_flex_self_energy();

// FFT-accelerated momentum-convolution variant of the same calculation, at
// the same Temperature. Requires k_mesh == q_mesh (see self_energy.cpp for
// why) and introduces a mesh-adaptive Lorentzian broadening eta instead of
// the exact analytic resolution - an approximation, but O(Nw^2 * Nk log Nk)
// instead of O(Nq*Nk*Nw). Only faster than call_flex_self_energy() when Nw
// is small relative to Nk (see self_energy.cpp for the crossover estimate).
void call_flex_self_energy_fft();
