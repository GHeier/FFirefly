# Just an example, not real
using BZIntegral.BZInt2D  # or BZInt3D

# Build meshes on your k-grid
Ek    = [ε(k)   for k in kgrid]
Ekpq  = [ε(k+q) for k in kgrid]
Zk    = [Z(k)   for k in kgrid]
Zkpq  = [Z(k+q) for k in kgrid]

# Denominator mesh: D(k) = ω - (ε(k+q) - ε(k)) + iη
Dk = (ω + 1im*η) .- Ekpq .+ Ek

# Smooth numerator F(k) = Z(k)*Z(k+q)
Fmesh = Zk .* Zkpq

# Term 1: Θ(eF - ε(k+q)) / D(k), weighted by F(k)
W1 = Quad2DRuleΘ𝔇(Ekpq, eF, Dk, iter)

# Term 2: Θ(eF - ε(k)) / D(k), weighted by F(k)
W2 = Quad2DRuleΘ𝔇(Ek, eF, Dk, iter)

χ = sum(Fmesh .* (W1 .- W2)) * vol
