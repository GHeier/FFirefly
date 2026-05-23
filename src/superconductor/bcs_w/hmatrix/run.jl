using Firefly
cfg = Firefly.Config

include("hmatrix_helper.jl")
using .HMatrixHelper

using HMatrices
using KrylovKit
using LinearAlgebra
using Printf
using StaticArrays

# Load relevant variables from the configuration
outdir = cfg.outdir
prefix = cfg.prefix
filetype = cfg.filetype

dim = cfg.dimension
Z = cfg.qp_weight
T = cfg.Temperature

mu = cfg.fermi_energy
n = cfg.num_electrons
if cfg.mu_from_n
    println("Initial mu = $mu")
    En = Firefly.Field_R(outdir * prefix * "_E_vs_n.h5")
    mu = En(n)
    println("Shifted mu = $mu")
end

# Must be even, code doesn't handle the pole well
w_pts = cfg.w_pts 
w_pts += (w_pts % 2) 

wc = cfg.cutoff_energy


function load_surface()
    eps_func = (k) -> Firefly.epsilon(1, [Float64(k.x), Float64(k.y), Float64(k.z)])

    # Create Surface at the Fermi energy (mu)
    surf = Firefly.Surface(eps_func, Float32(mu))
    kpoints, areas = get_faces_and_areas(surf)
    println("Found $(length(kpoints)) k-points on Fermi surface")
    println("Total area: $(sum(areas))")

    return kpoints, areas
end

function get_surface_data()
    # Parameters
    println("\nGenerating k-points from Fermi surface at μ = $mu")
    H = Firefly.Hamiltonian()

    w_points = Vector{Float32}(undef, w_pts)
    if w_pts == 1
        # Single frequency point: use w=0 (static pairing)
        w_points[1] = 0.0
    else
        # Multiple frequency points: distribute from 0 to wc
        for i in 1:w_pts
            w_points[i] = -wc + (2 * wc) * (i - 1) / (w_pts - 1)
            #w_points[i] = wc * (i - 1) / (w_pts - 1)
        end
    end


    kpoints, areas = load_surface()
    n = length(kpoints)
    println("Total k-points: $n")

    # Load Fermi velocities at k-points
    println("\nComputing Fermi velocities...")
    get_fermi_velocity = Firefly.Imports.get_fermi_velocity
    # Convert Float32 k-points to Float64 for get_fermi_velocity
    kpoints_f64 = [Float64.(kp) for kp in kpoints]
    velocities = get_fermi_velocity(H, kpoints_f64)
    v_norms = [LinearAlgebra.norm(velocities[i, 1, :]) for i in 1:n]
    @printf("Fermi Velocity Min, Max, Ave: %.6f, %.6f, %.6f\n", minimum(v_norms), maximum(v_norms), sum(v_norms) / n)

    # Compute density of states weights
    dos_weights = areas ./ v_norms ./ (2 * π)^dim
    @printf("DOS Min, Max, Ave, Total: %.6f, %.6f, %.6f, %.6f\n", minimum(dos_weights), maximum(dos_weights), sum(dos_weights) / n, sum(dos_weights))

    return kpoints, dos_weights, w_points
end

function f(w)
    return 1 / (exp(w / T) + 1)
end

function make_vertex_kernel(V, weights, frequencies, with_w)
    # Return kernel function V(k-k') for k-k' only (no frequency dependence)
    return function(k1, k2, i, j)
        weight1 = weights[i]
        weight2 = weights[j]
        k1 = Float64.(k1)
        k2 = Float64.(k2)

        dkm = k1 - k2
        dkp = k1 + k2
        V_val = real(V(dkm) + V(dkp)) / 2
        if with_w                                                                               
            dV = 0
            for w in frequencies
                if (w) < 1e-4
                    continue
                end
                dV += real( (V(dkm,w) + V(dkp,w)) / 2 - V_val) / w
            end
            V_val = dV
        end
        #V_val = 1.0

        return -(weight1 * weight2)^0.5 * V_val
    end
end



function create_hmatrices(V, dos_weights, kpoints, frequencies)
    # Create kernel function V(k-k') for k-k' only
    kernel = make_vertex_kernel(V, dos_weights, frequencies, false)
    kernel_ptb = make_vertex_kernel(V, dos_weights, frequencies, true)

    # Build HMatrix for k-k' points only (no frequency dependence)
    println("\n" * "="^60)
    println("\nBuilding HMatrix for V(k-k')...")
    t_hmat = @elapsed H_matrix = HMatrixHelper.build_hmatrix(kpoints, kernel; atol=1e-4, rank=30)
    @printf("HMatrix built in %.3f seconds\n", t_hmat)
    @printf("Compression ratio: %.2f\n", (compression_ratio(H_matrix)))
    println("\n" * "="^60)
    println("\nBuilding HMatrix for V(k-k') Perturbation...")
    t_hmat_ptb = @elapsed H_matrix_ptb = HMatrixHelper.build_hmatrix(kpoints, kernel_ptb; atol=1e-4, rank=30)
    @printf("HMatrix built in %.3f seconds\n", t_hmat_ptb)
    @printf("Compression ratio: %.2f\n", (compression_ratio(H_matrix_ptb)))

    # Probe a few raw kernel values to verify non-zero perturbation
    n_k = length(kpoints)
    println("\nDEBUG: Sample kernel values (i,j) for i!=j:")
    for (i,j) in [(1,2),(1,3),(n_k÷2, n_k÷2+1)]
        v0 = kernel(kpoints[i], kpoints[j], i, j)
        vp = kernel_ptb(kpoints[i], kpoints[j], i, j)
        @printf("  (%d,%d): K0=%.6e  Kptb=%.6e  ratio=%.4f\n", i, j, v0, vp, abs(vp)/(abs(v0)+1e-30))
    end

    return H_matrix, H_matrix_ptb
end

function lanczos_solve(H_matrix, kpoints)
    # Matrix dimensions: n_k k-points only (frequency factor applied as scalar)
    n_k = length(kpoints)
    @printf("Matrix size: %d × %d\n", n_k, n_k)

    hmv = (v) -> begin
        result = similar(v)
        mul!(result, H_matrix, v)
        return result
    end

    num_eigs = min(10, n_k)  # Request up to 10 eigenvalues
    t_lanczos_hmat = @elapsed begin # Perform Lanczos
        vals_hmat, vecs_hmat, info_hmat = eigsolve(hmv, n_k, num_eigs, :LR; issymmetric=false, krylovdim=30, maxiter=200, tol=1e-8)
    end

    return vals_hmat, vecs_hmat, info_hmat, t_lanczos_hmat
end

function perturb_solve(H, eigs, vecs)
    hmv_eig = (v) -> begin
        result = similar(v)
        mul!(result, H, v)
        return sum(result .* v)
    end

    hmv_vec = (v1, v2) -> begin
        result = similar(v1)
        mul!(result, H, v1)
        l = sum(result .* v2)
        #println("Project val: ", l)
        return  l .* v2 
    end

    # --- Projection diagnostics for dominant mode (i=1) ---
    v0 = real.(vecs[1])
    Hv0 = similar(v0)
    mul!(Hv0, H, v0)  # H_ptb * v0
    println("\n=== Projection of H_ptb * v0 onto eigenbasis of H_0 ===")
    @printf("  Diagonal <v0|H_ptb|v0>     = %.6e  (this is λ1)\n", LinearAlgebra.dot(v0, Hv0))
    println("  Off-diagonal <vj|H_ptb|v0>:")
    offdiag_norm_sq = 0.0
    for j in 2:length(vecs)
        vj = real.(vecs[j])
        c = LinearAlgebra.dot(vj, Hv0)
        eig_gap = real(eigs[1] - eigs[j])
        @printf("    j=%2d: <vj|H_ptb|v0>=%.4e  eig_gap=%.4e  correction_coeff=%.4e\n",
            j, c, eig_gap, c / eig_gap)
        offdiag_norm_sq += c^2
    end
    @printf("  |off-diagonal part of H_ptb*v0| = %.6e\n", sqrt(offdiag_norm_sq))
    @printf("  |H_ptb*v0| total               = %.6e\n", LinearAlgebra.norm(Hv0))
    @printf("  Fraction captured by 10 eigvecs = %.4f\n",
        sqrt(LinearAlgebra.dot(v0,Hv0)^2 + offdiag_norm_sq) / LinearAlgebra.norm(Hv0))
    println("======================================================")

    ptb_eigs = zeros(Float64, length(eigs))
    ptb_vecs = [zeros(Float64, length(v)) for v in vecs]
    for i in 1:length(vecs)
        ptb_eigs[i] = real(hmv_eig(vecs[i]))
        ptb_vecs[i] .= 0.0
        for j in 1:length(vecs)
            i == j && continue
            eig_gap = eigs[i] - eigs[j]
            proj = hmv_vec(vecs[i], vecs[j])
            correction = real.(proj / eig_gap)
            ptb_vecs[i] += correction
        end
    end
    @printf("DEBUG: |ptb_vecs[1]|=%.6e, |vecs[1]|=%.6e\n", LinearAlgebra.norm(ptb_vecs[1]), LinearAlgebra.norm(vecs[1]))
    @printf("DEBUG: overlap(ptb_vecs[1], vecs[1])=%.6e\n", LinearAlgebra.dot(ptb_vecs[1], real.(vecs[1])))

    return ptb_eigs, ptb_vecs
end

function find_top_eig(vals_hmat, vecs_hmat, vals_ptb, vecs_ptb)
    l0 = real.(vals_hmat)
    l1 = real.(vals_ptb)

    eff_lambdas = zeros(size(vals_hmat))
    eff_vecs = [copy(v) for v in vecs_hmat]  # deep copy each eigenvector
    mu_star = 0.10
    @printf("\nTop eigenvalues:\n")
    for i in 1:length(eff_lambdas)
        v0 = real.(vecs_hmat[i])  # snapshot BEFORE mutation
        eff_lambdas[i] = l0[i] * (1 - mu_star) / (1 / Z + l1[i])
        eff_vecs[i] .+= vecs_ptb[i]
        vptb = vecs_ptb[i]
        eff_v = eff_vecs[i]
        @printf("%d)  λ0=%.8f λ1=%.8f λe=%.8f  |ptb_vec|=%.4e  cos(v0,eff)=%.6f\n",
            i, l0[i], l1[i], eff_lambdas[i], LinearAlgebra.norm(vptb),
            LinearAlgebra.dot(v0, real.(eff_v)) / (LinearAlgebra.norm(v0) * LinearAlgebra.norm(real.(eff_v)) + 1e-30))
    end
    println("Quasiparticle Weight Z = ", Z)
    println("mu* = ", mu_star)
    lz = 1/Z - 1
    @printf("λz = %.10f\n", lz)
    perm = sortperm(eff_lambdas)
    l_sort = eff_lambdas[perm]
    v_sort = eff_vecs[perm]
    l0_sort = l0[perm]
    l1_sort = l1[perm]

    if l_sort[end] <= 0
        println("No positive eigenvalues found :(")
        exit(1)
    end

    return l_sort, v_sort, l0_sort, l1_sort
end

function un_hermitize_gap(gap, weights)
    # Un-hermitize by dividing by sqrt(weights)
    # Note: weights include (2π)^dim factor for hermitization, but we need to
    # remove it here to match C++ scaling which uses sqrt(area/v_norm) directly
    result = copy(gap)
    factor = (2 * π)^(dim / 2)  # dim is global
    for i in 1:length(result)
        result[i] /= weights[i]^0.5 * factor
    end
    return Float64.(real.(result))
end

function save!(vals, vecs, kpoints, weights)
    n_k = length(kpoints)
    k_dim = length(kpoints[1])
    # Prepare k-points matrix [n_k × k_dim]
    points_matrix = Matrix{Float64}(undef, n_k, k_dim)
    for i_k in 1:n_k
        for d in 1:k_dim
            points_matrix[i_k, d] = Float64(kpoints[i_k][d])
        end
    end

    save_data! = Firefly.Imports.save_data!
    for i in 1:length(vals)
        filename = outdir * prefix * "_gap_$i." * filetype
        @printf("Saving eigenvector to %s\n", filename)
        data_array = un_hermitize_gap(vecs[i], weights)
        save_data!(filename, data_array, points=points_matrix)
    end
    # KrylovKit's eigsolve with :LR returns eigenvalues sorted largest-first,
    # so vecs[1] is the dominant eigenvector (largest eigenvalue), not vecs[end]
    filename = outdir * prefix * "_gap.h5"
    save_data!(filename, un_hermitize_gap(vecs[end], weights), points=points_matrix)
    @printf("Saved dominant eigenvector to %s\n", filename)
end

function run()
    # Main function call goes here
    println("="^60)
    println("Initializing HMatrix + Arnoldi Eigenvalue Solver (k-k' only)")
    println("="^60)

    kpoints, dos_weights, w_points = get_surface_data()

    # Load Vertex
    println("\nLoading Vertex...")
    V = Firefly.Field_R(outdir * prefix * "_vertex_singlet.h5")
    println("Vertex loaded.")

    H_matrix, H_matrix_ptb = create_hmatrices(V, dos_weights, kpoints, w_points)
    println("Created HMatrix")
    vals_hmat, vecs_hmat, info_hmat, t_lanczos_hmat = lanczos_solve(H_matrix, kpoints)
    #println(sort(real.(vals_hmat), rev=true))
    println("Performing BCS+w perturbation")
    ptb_eigs, ptb_vecs = perturb_solve(H_matrix_ptb, vals_hmat, vecs_hmat)

    @printf("Time: %.3f seconds\n", t_lanczos_hmat)
    @printf("Number of eigenvalues computed: %d\n", length(vals_hmat))
    @printf("Converged: %s (iterations: %d)\n", info_hmat.converged > 0, info_hmat.numiter)

    l, v, l0, l1 = find_top_eig(vals_hmat, vecs_hmat, ptb_eigs, ptb_vecs)
    println("Max Eig: ", l[end])

    fT = log(1.134 * wc / T)
    println("Max f(T)*Eig: ", fT * l[end])
    Tc = 1.134 * wc * exp(-1/l[end])
    println("Tc (eV): ", Tc)
    println("Tc (K): ", Tc * 11604.5)

    save!(l, v, kpoints, dos_weights)

    return l[end]
end

if abspath(PROGRAM_FILE) == @__FILE__ # Runs on file execution
    run()
end



