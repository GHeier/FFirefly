#!/usr/bin/env julia

using Firefly

# Load configuration
Firefly.load_config!("/home/g/Research/Materials/Tight_Binding/test/sample.cfg")

println("="^60)
println("Testing Surface with get_faces_and_areas in Julia")
println("="^60)

# Create epsilon function
eps_func = (k) -> Firefly.epsilon(1, [Float64(k.x), Float64(k.y), Float64(k.z)])

# Create Surface at Fermi energy
mu = -1.5
println("\nCreating surface at μ = $mu")
surf = Firefly.Surface(eps_func, Float32(mu))

# Test old method
println("\nTesting get_faces (old method):")
faces = Firefly.get_faces(surf)
println("  Number of faces: $(length(faces))")
println("  First face k-point: $(faces[1])")
println("  Dimension of first face: $(length(faces[1]))")

# Test new method
println("\nTesting get_faces_and_areas (new method):")
kpoints, areas = Firefly.get_faces_and_areas(surf)
println("  Number of k-points: $(length(kpoints))")
println("  Number of areas: $(length(areas))")
println("  First k-point: $(kpoints[1])")
println("  First area: $(areas[1])")

# Verify they match
println("\nVerification:")
println("  K-points match: $(kpoints == faces)")
println("  Total area: $(sum(areas))")
println("  Min area: $(minimum(areas))")
println("  Max area: $(maximum(areas))")
println("  Mean area: $(sum(areas)/length(areas))")

println("\n" * "="^60)
println("Julia test completed successfully!")
println("="^60)
