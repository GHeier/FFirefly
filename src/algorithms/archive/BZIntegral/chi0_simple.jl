using BZIntegral
using BZIntegral.BZInt2D
using LinearAlgebra
using Printf

"""
Bare susceptibility χ₀(q) calculation for a 2D tight binding model.

This script uses the BZIntegral.jl package to accurately compute the
Lindhard function (bare electronic susceptibility) at zero temperature
and zero frequency for a square lattice tight binding model.

Model:
  ε(k) = -2t[cos(kx) + cos(ky)] - μ
  Parameters: t = 1.0, μ = -1.0

Method:
  Uses the recursive hybrid tetrahedron method from BZIntegral.jl
  to handle the singular 1/(εₖ₊q - εₖ) denominator.
"""

function tight_binding_energy_2d(kx, ky, t=1.0, mu=-1.0)
    """2D tight binding dispersion"""
    return -2.0 * t * (cos(kx) + cos(ky)) - mu
end

function calculate_chi0_simple(qx, qy, nk=120, t=1.0, mu=-1.0, iter=2)
    """
    Calculate the bare susceptibility χ₀(q) at zero temperature and frequency.

    The susceptibility is:
        χ₀(q,ω=0) = -∫ d²k/(2π)² [f(εₖ) - f(εₖ₊q)] / (εₖ₊q - εₖ)

    At T=0: f(ε) = θ(μ - ε), where θ is the Heaviside step function.

    Arguments:
        qx, qy: Momentum transfer components
        nk: Number of k-points per dimension (default: 120)
        t: Hopping parameter (default: 1.0)
        mu: Chemical potential (default: -1.0)
        iter: BZIntegral recursion depth (default: 2)

    Returns:
        χ₀(q): Bare susceptibility at the given momentum transfer
    """

    # Create meshes for k-space
    # For tight binding, we integrate over [-π, π] × [-π, π]
    kx = range(-π, π, length=nk)
    ky = range(-π, π, length=nk)

    # Create 2D meshes
    KX = repeat(kx, 1, nk)
    KY = repeat(ky', nk, 1)

    # Calculate dispersions
    Ek = zeros(nk, nk)
    Ekplusq = zeros(nk, nk)

    for i in 1:nk
        for j in 1:nk
            Ek[i,j] = tight_binding_energy_2d(KX[i,j], KY[i,j], t, mu)
            # k+q with periodic BCs
            kxq = mod(KX[i,j] + qx + π, 2π) - π
            kyq = mod(KY[i,j] + qy + π, 2π) - π
            Ekplusq[i,j] = tight_binding_energy_2d(kxq, kyq, t, mu)
        end
    end

    # Fermi level (energies are measured from μ, so eF = 0)
    eF = 0.0
    ω = 0.0

    # Denominator: (ω + εₖ - εₖ₊q)
    # Try without the 2π factor first to see if it's a normalization issue
    Dk = (ω .+ Ek .- Ekplusq)

    # Convert to OBC (open boundary conditions) as required by BZIntegral
    Ek_obc = PBC2OBC_2D(Ek)
    Ekplusq_obc = PBC2OBC_2D(Ekplusq)
    Dk_obc = PBC2OBC_2D(Dk)

    # Calculate using Quad2DRuleΘ𝔇
    # This computes: ∫ θ(eF - E(k)) / D(k)
    Wmesh_obc = Quad2DRuleΘ𝔇(Ek_obc, eF, Dk_obc, iter) - Quad2DRuleΘ𝔇(Ekplusq_obc, eF, Dk_obc, iter)

    # Convert back to PBC
    Wmesh = OBC2PBC_2D(Wmesh_obc)

    # Integrate: sum over mesh points
    # The BZIntegral tetrahedron method returns weights that when summed
    # give the integral over the unit volume. We need to scale by BZ volume.
    chi0 = sum(Wmesh)

    return chi0
end

function test_chi0_path()
    """Test χ₀ along (0,0) → (π,π) path"""
    println("=" ^ 70)
    println("2D Tight Binding Bare Susceptibility χ₀(q)")
    println("=" ^ 70)
    println("Parameters: t = 1.0, μ = -1.0")
    println("Dispersion: ε(k) = -2t[cos(kx) + cos(ky)] - μ")
    println()

    nk = 120
    nq = 11

    println(@sprintf("k-mesh: %d × %d", nk, nk))
    println(@sprintf("q-path: (0,0) → (π,π) with %d points", nq))
    println()
    println("-" ^ 70)
    println(@sprintf("%10s %15s %15s %15s", "Point", "qx", "qy", "χ₀(q)"))
    println("-" ^ 70)

    results = []

    for i in 0:(nq-1)
        t_param = i / (nq - 1)
        qx = t_param * π
        qy = t_param * π

        # Skip q=0 as it diverges
        if i == 0
            println(@sprintf("%10d %15.4f %15.4f %15s", i, qx, qy, "(divergent)"))
            continue
        end

        chi0 = calculate_chi0_simple(qx, qy, nk)

        push!(results, (qx, qy, chi0))
        println(@sprintf("%10d %15.4f %15.4f %15.6f", i, qx, qy, chi0))
    end

    println("-" ^ 70)
    println()

    # Check results
    chi_values = [r[3] for r in results]
    chi_min, chi_max = extrema(chi_values)

    println("Summary:")
    println(@sprintf("  Min χ₀: %.6f", chi_min))
    println(@sprintf("  Max χ₀: %.6f", chi_max))

    # Check magnitude (susceptibility is negative by convention in some definitions)
    abs_min = abs(chi_min)
    abs_max = abs(chi_max)

    if abs_min >= 0.13 && abs_max <= 0.22
        println("  ✓ Results in expected range [0.13, 0.22]")
    elseif abs_min >= 0.10 && abs_max <= 0.25
        println("  ⚠ Results close to expected range:")
        println("    |χ₀| ∈ [$(round(abs_min, digits=4)), $(round(abs_max, digits=4))]")
        println("    Expected: [0.13, 0.22]")
    else
        println("  Note: |χ₀| ∈ [$(round(abs_min, digits=4)), $(round(abs_max, digits=4))]")
        println("  (Expected: 0.13 to 0.22)")
    end
    println()

    return results
end

println("\nStarting calculation...\n")
try
    results = test_chi0_path()
    println("Calculation completed!")
catch e
    println("Error: $e")
    rethrow(e)
end
