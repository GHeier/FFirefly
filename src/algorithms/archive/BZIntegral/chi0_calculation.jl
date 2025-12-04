using BZIntegral
using BZIntegral.BZInt2D
using LinearAlgebra
using Printf

"""
Calculate the bare susceptibility χ₀(q) for a 2D tight binding model
using the BZIntegral.jl package.

Model: ε(k) = -2t(cos(kx) + cos(ky)) - μ
Parameters: t = 1.0, μ = -1.0
"""

function tight_binding_energy_2d(kx, ky, t=1.0, mu=-1.0)
    """2D tight binding dispersion"""
    return -2.0 * t * (cos(kx) + cos(ky)) - mu
end

function calculate_chi0_2d(qx, qy, nk=100, t=1.0, mu=-1.0, T=1e-6, iter=2)
    """
    Calculate bare susceptibility χ₀(q) at a given momentum transfer q=(qx,qy)

    For zero temperature and ω=0, the Lindhard function is:
    χ₀(q) = -∑ₖ [f(εₖ) - f(εₖ₊q)] / (εₖ₊q - εₖ)

    Following the BZIntegral.jl Lindhard2D.jl example:
    Wmesh = Quad2DRuleΘ𝔇(Ek,eF,Dk) - Quad2DRuleΘ𝔇(Ekplusq,eF,Dk)
    where Dk = (ω + Ek - Ekplusq) * 2π
    """

    # Create k-mesh in the first Brillouin zone [-π, π]
    kx_grid = range(-π, π, length=nk)
    ky_grid = range(-π, π, length=nk)

    # Calculate energy meshes
    Emesh_k = zeros(nk, nk)
    Emesh_kpq = zeros(nk, nk)

    for (i, kx) in enumerate(kx_grid)
        for (j, ky) in enumerate(ky_grid)
            ek = tight_binding_energy_2d(kx, ky, t, mu)
            # k+q with periodic boundary conditions
            kx_q = mod(kx + qx + π, 2π) - π
            ky_q = mod(ky + qy + π, 2π) - π
            ekpq = tight_binding_energy_2d(kx_q, ky_q, t, mu)

            Emesh_k[i, j] = ek
            Emesh_kpq[i, j] = ekpq
        end
    end

    # Fermi energy (we measure from chemical potential, so eF = 0)
    eF = 0.0
    ω = 0.0  # Zero frequency

    # Calculate Dmesh = (ω + εₖ - εₖ₊q) * 2π
    # Note the sign: Dk = (ω + Ek - Ekplusq) from the example
    Dmesh = (ω .+ Emesh_k .- Emesh_kpq) .* (2π)

    # Convert to OBC for BZIntegral
    Emesh_k_obc = PBC2OBC_2D(Emesh_k)
    Emesh_kpq_obc = PBC2OBC_2D(Emesh_kpq)
    Dmesh_obc = PBC2OBC_2D(Dmesh)

    # Calculate the two terms
    Wmesh1_obc = Quad2DRuleΘ𝔇(Emesh_k_obc, eF, Dmesh_obc, iter)
    Wmesh2_obc = Quad2DRuleΘ𝔇(Emesh_kpq_obc, eF, Dmesh_obc, iter)

    # The Wmesh values need to be converted back
    Wmesh1 = OBC2PBC_2D(Wmesh1_obc)
    Wmesh2 = OBC2PBC_2D(Wmesh2_obc)

    # Sum over BZ
    # The BZIntegral example uses: out = sum(Wmesh) * vol
    # where vol is the actual volume of the integration region
    vol = (2π)^2 / (nk^2)
    chi0 = sum(Wmesh1 .- Wmesh2) * vol

    return chi0
end

function test_chi0_along_path()
    """
    Test χ₀ along the path from (0,0) to (π,π)
    Expected values: between 0.13 and 0.22
    """
    println("=" ^ 70)
    println("Bare Susceptibility χ₀(q) for 2D Tight Binding Model")
    println("=" ^ 70)
    println("Parameters: t = 1.0, μ = -1.0, T → 0")
    println("Model: ε(k) = -2t(cos(kx) + cos(ky)) - μ")
    println()

    # Number of q-points along the path
    nq = 11
    nk = 120  # k-mesh resolution

    println(@sprintf("k-mesh: %d × %d points", nk, nk))
    println(@sprintf("q-path: (0,0) → (π,π) with %d points", nq))
    println()
    println("-" ^ 70)
    println(@sprintf("%10s %15s %15s %15s", "Point", "qx", "qy", "χ₀(q)"))
    println("-" ^ 70)

    results = []

    for i in 0:(nq-1)
        # Path from (0,0) to (π,π)
        t_param = i / (nq - 1)
        qx = t_param * π
        qy = t_param * π

        # Calculate χ₀
        chi0 = calculate_chi0_2d(qx, qy, nk)

        push!(results, (qx, qy, chi0))
        println(@sprintf("%10d %15.4f %15.4f %15.6f", i, qx, qy, chi0))
    end

    println("-" ^ 70)
    println()

    # Check if results are in expected range
    chi_values = [r[3] for r in results]
    chi_min, chi_max = extrema(chi_values)

    println("Summary:")
    println(@sprintf("  Min χ₀: %.6f", chi_min))
    println(@sprintf("  Max χ₀: %.6f", chi_max))

    if chi_min >= 0.13 && chi_max <= 0.22
        println("  ✓ Results are in expected range [0.13, 0.22]")
    elseif chi_min >= 0.10 && chi_max <= 0.30
        println("  ⚠ Results are close to expected range [0.13, 0.22]")
        println("    (obtained: [$(round(chi_min, digits=3)), $(round(chi_max, digits=3))])")
    else
        println("  ✗ Results are outside expected range [0.13, 0.22]")
        println("    (obtained: [$(round(chi_min, digits=3)), $(round(chi_max, digits=3))])")
    end
    println()

    return results
end

# Run the test
println("\nStarting calculation...")
println("This may take a few moments...\n")

try
    results = test_chi0_along_path()
    println("Calculation completed successfully!")
catch e
    println("Error occurred during calculation:")
    println(e)
    rethrow(e)
end
