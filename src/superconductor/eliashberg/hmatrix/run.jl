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
T = cfg.Temperature
Z = cfg.qp_weight

mu = cfg.fermi_energy
n = cfg.num_electrons
mu_from_n = cfg.mu_from_n
#
# Must be even, code doesn't handle the pole well
w_pts = cfg.w_pts 
w_pts += (w_pts % 2) 

wc = cfg.cutoff_energy

debug = false

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
    println("")


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

function make_vertex_kernel(V, weights, frequencies::Vector{Float32})
    dw = 2 * wc / w_pts
    # Return kernel function V(k-k', w+w')
    # Signature: kernel(k1, k2, w1, w2, i, j, iw, jw)
    return function(k1, k2, w1, w2, i, j, iw, jw)
        # Convert to Float64 for Vertex call
        k1 = Float64.(k1)
        k2 = Float64.(k2)
        w1 = Float64(w1)
        w2 = Float64(w2)

        dk = k1 - k2
        V_val = 1.0
        if !debug
            V_val = real(V(dk, w1 + w2))
        end
        weight = weights[j]

        return dw * weight * V_val
    end
end



function create_hmatrices(V, dos_weights, kpoints, w_points)
    # Create two kernel functions: V(k-k', w-w') and V(k-k', w+w')
    kernel = make_vertex_kernel(V, dos_weights, w_points)
    total_N = size(dos_weights,1) * size(w_points,1)
    println("Total N = $total_N")

    # Build two HMatrices with frequency dependence
    println("\n" * "="^60)
    println("\nBuilding HMatrix for V(k-k', w+w')...")
    t_hmat= @elapsed H_matrix = HMatrixHelper.build_hmatrix(kpoints, w_points, kernel; atol=1e-4, rank=30)
    @printf("HMatrices built in %.3f seconds\n", t_hmat)
    @printf("Compression ratio: %.2f\n", (compression_ratio(H_matrix)))

    return H_matrix
end

function lanczos_solve(H_matrix, kpoints, w_points)
    # Matrix dimensions: n_k k-points × n_w frequencies
    n_k = length(kpoints)
    n_w = length(w_points)
    n_total = n_k * n_w
    @printf("Matrix size: %d × %d\n", n_total, n_total)

    hmv = (v) -> begin
        newv = reshape(v, n_k, n_w)  # [k, w]
        fv = similar(newv)

        for i_w in 1:n_w
            w = Float64(w_points[i_w])
            if abs(w) < 1e-8
                fv[:, i_w] = newv[:, i_w] .* (-1 / (8 * T))
            else
                fv[:, i_w] = newv[:, i_w] .* (-f(w) / (w))
            end
        end

        result = similar(v)
        mul!(result, H_matrix, reshape(fv, n_total))
        #result .= 1 / n_total

        return result
    end

    num_eigs = min(10, n_total)  # Request up to 10 eigenvalues
    t_lanczos_hmat = @elapsed begin # Perform Lanczos
        vals_hmat, vecs_hmat, info_hmat = eigsolve(hmv, n_total, num_eigs, :LR; issymmetric=false, krylovdim=30, maxiter=200, tol=1e-8)
    end

    return vals_hmat, vecs_hmat, info_hmat, t_lanczos_hmat
end

function find_top_eig(vals_hmat, vecs_hmat, info_hmat, t_lanczos_hmat)
    real_vals = real.(vals_hmat)
    positive_vals = filter(x -> x > 0, real_vals)

    if isempty(positive_vals)
        @printf("Warning: No positive eigenvalues found!\n")
        @printf("All %d eigenvalues:\n", length(real_vals))
        for (i, val) in enumerate(real_vals)
            @printf("  λ_%d = %.10f\n", i, val)
        end
        largest_positive = NaN
        idx_largest_positive = 0
    else
        largest_positive = maximum(positive_vals)
        # Find the index of the largest positive eigenvalue
        idx_largest_positive = findfirst(x -> x == largest_positive, real_vals)
    end

    @printf("Time: %.3f seconds\n", t_lanczos_hmat)
    @printf("Number of eigenvalues computed: %d\n", length(vals_hmat))
    @printf("Converged: %s (iterations: %d)\n", info_hmat.converged > 0, info_hmat.numiter)

    # Show top eigenvalues
    @printf("\nTop eigenvalues:\n")
    for i in 1:min(5, length(real_vals))
        @printf("  λ_%d = %.10f\n", i, real_vals[i])
    end

    @printf("\nLargest eigenvalue (by magnitude): %.10f\n", real_vals[1])
    @printf("Largest positive eigenvalue: %.10f\n", largest_positive)
    return idx_largest_positive, largest_positive
end

# Save the eigenvector corresponding to the largest positive eigenvalue
function save!(vals, vecs, idx, kpoints, w_points)
    if idx > 0
        eigenvec = vecs[idx]
        n_k = length(kpoints)
        n_w = length(w_points)

        eigenvec_2d = reshape(eigenvec, n_k, n_w)

        # Save gap function with k-points from Fermi surface using save_data!
        filename = outdir * prefix * "_gap." * filetype
        @printf("\nSaving eigenvector to %s\n", filename)

        # Determine k-space dimension
        k_dim = length(kpoints[1])

        # Prepare k-points matrix [n_k × k_dim]
        points_matrix = Matrix{Float64}(undef, n_k, k_dim)
        for i_k in 1:n_k
            for d in 1:k_dim
                points_matrix[i_k, d] = Float64(kpoints[i_k][d])
            end
        end

        data_array = Float64.(real.(eigenvec_2d))
        save_data! = Firefly.Imports.save_data!
        if w_pts > 1
            save_data!(filename, data_array, points=points_matrix, w_points=w_points)
        else
            save_data!(filename, data_array, points=points_matrix)
        end

        @printf("Eigenvector saved successfully.\n")
    else
        @printf("\nNo positive eigenvalue found - eigenvector not saved.\n")
    end
end

function run()
    global mu
    global mu_from_n
    if mu_from_n
        println("Initial mu = $mu")
        En = Firefly.Field_R(outdir * prefix * "_E_vs_n.h5")
        mu = En(n)
        println("Shifted mu = $mu")
    end

    if debug
        println("Running in DEBUG mode")
    end
    # Main function call goes here
    println("="^60)
    println("Initializing HMatrix + Arnoldi Eigenvalue Solver")
    println("="^60)

    kpoints, dos_weights, w_points = get_surface_data()

    println("Quasiparticle Weight: ", Z)

    if !debug
        # Load Vertex
        println("\nLoading Vertex...")
        V = Firefly.Field_C(outdir * prefix * "_vertex_singlet.h5")
        println("Vertex loaded.")
    else
        V = 1
    end

    H_matrix = create_hmatrices(V, dos_weights, kpoints, w_points)
    println("Created HMatrix")
    vals_hmat, vecs_hmat, info_hmat, t_lanczos_hmat = lanczos_solve(H_matrix, kpoints, w_points)
    if !debug
        vals_hmat .*= Z
    end

    fT = log(1.134 * wc / T)
    eig_est = fT * sum(dos_weights)
    println("Eigenvalue estimate: ", eig_est)

    idx_largest_positive, largest_positive = find_top_eig(vals_hmat, vecs_hmat, info_hmat, t_lanczos_hmat)
    save!(vals_hmat, vecs_hmat, idx_largest_positive, kpoints, w_points)

    return largest_positive, eig_est
end

if abspath(PROGRAM_FILE) == @__FILE__ # Runs on file execution
    run()
end



