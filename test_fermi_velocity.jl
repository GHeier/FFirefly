#!/usr/bin/env julia

using Firefly
using LinearAlgebra: norm
using Printf
# Use Imports module directly for new functions
const get_fermi_velocity = Firefly.Imports.get_fermi_velocity

println("="^60)
println("Testing Hamiltonian.get_fermi_velocity in Julia")
println("="^60)

# Load configuration
Firefly.load_config!("/home/g/Research/Materials/Tight_Binding/test/sample.cfg")

# Create Hamiltonian
println("\nInitializing Hamiltonian...")
H = Firefly.Hamiltonian()
println("Hamiltonian loaded from file: ", Firefly.file_found(H))

# Test 1: Single k-point
println("\n" * "="^60)
println("Test 1: Single k-point")
println("="^60)

k = [0.1, 0.2, 0.0]
println("k-point: $k")

# Get bands
bands = Firefly.get_bands(H, k)
println("Number of bands: $(length(bands))")
println("Bands: $bands")

# Get Fermi velocity
vels = get_fermi_velocity(H, k)
println("\nFermi velocities (size: $(size(vels))):")
for i in 1:size(vels, 1)
    v = vels[i, :]
    v_norm = norm(v)
    @printf("  Band %d: v = [%.6f, %.6f, %.6f], |v| = %.6f\n", i, v[1], v[2], v[3], v_norm)
end

# Test 2: Multiple k-points
println("\n" * "="^60)
println("Test 2: Multiple k-points")
println("="^60)

k_points = [
    [0.0, 0.0, 0.0],
    [0.5, 0.0, 0.0],
    [0.5, 0.5, 0.0],
    [0.0, 0.5, 0.0]
]
println("Number of k-points: $(length(k_points))")

# Get velocities for all k-points
vels_list = get_fermi_velocity(H, k_points)
println("Velocities array size: $(size(vels_list))")

for (i, k) in enumerate(k_points)
    println("\nk[$i] = $k")
    for n in 1:size(vels_list, 2)
        v = vels_list[i, n, :]
        v_norm = norm(v)
        @printf("  Band %d: v = [%.4f, %.4f, %.4f], |v| = %.4f\n", n, v[1], v[2], v[3], v_norm)
    end
end

# Test 3: Verify numerical derivative is correct
println("\n" * "="^60)
println("Test 3: Verify numerical derivative")
println("="^60)

k0 = [0.3, 0.4, 0.0]
dk = 0.001

println("Testing at k = $k0")
println("Using finite difference step dk = $dk")

# Get Fermi velocity from our function
v_computed = get_fermi_velocity(H, k0)

# Compute numerical derivative manually
E0 = Firefly.get_bands(H, k0)
nbands = length(E0)

kx_plus = [k0[1] + dk, k0[2], k0[3]]
ky_plus = [k0[1], k0[2] + dk, k0[3]]

Ex_plus = Firefly.get_bands(H, kx_plus)
Ey_plus = Firefly.get_bands(H, ky_plus)

v_manual = zeros(Float32, nbands, 3)
for n in 1:nbands
    v_manual[n, 1] = (Ex_plus[n] - E0[n]) / dk
    v_manual[n, 2] = (Ey_plus[n] - E0[n]) / dk
    v_manual[n, 3] = 0.0
end

println("\nComparison:")
for n in 1:nbands
    v_comp = v_computed[n, :]
    v_man = v_manual[n, :]
    diff = norm(v_comp - v_man)
    @printf("  Band %d:\n", n)
    @printf("    Computed: [%.6f, %.6f, %.6f]\n", v_comp[1], v_comp[2], v_comp[3])
    @printf("    Manual:   [%.6f, %.6f, %.6f]\n", v_man[1], v_man[2], v_man[3])
    @printf("    Difference: %.8f\n", diff)
end

println("\n" * "="^60)
println("Julia test completed successfully!")
println("="^60)
