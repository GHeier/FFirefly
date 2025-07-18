module Eliashberg
println("Started Julia")
using Firefly

t1 = time()
include("../objects/mesh.jl")
using .IRMesh

using Printf
using CUDA, FFTW
using MPI
#using NFFT
#using Plots
using Roots, Statistics
using SparseIR
import SparseIR: Statistics, value, valueim, MatsubaraSampling64F, TauSampling64
using Base.Threads
#using JLD
using LinearAlgebra, Printf, PyCall
np = pyimport("numpy")

t2 = time()
#println("Time to load: ", t2 - t1)

t3 = time()
#println("Time to load ffirefly: ", t3 - t2)
cfg = Firefly.Config

prefix = cfg.prefix
outdir = cfg.outdir

np_BZ = np.array(cfg.brillouin_zone)
projections = cfg.projections

nx, ny, nz = cfg.k_mesh
dim = cfg.dimension
#println("Dimension: ", dim)
if dim == 2
    nz = 1
end
nk = nx * ny * nz
nw = cfg.w_pts
U = cfg.onsite_U
BZ = cfg.brillouin_zone
mu = cfg.fermi_energy
t4 = time()
#println("Time to load config: ", t4 - t3)

const beta = 1 / cfg.Temperature
println("Beta: ", beta)
const pi = π

wD = cfg.cutoff_energy

function to_IBZ(k)
    tolerance = 1e-10  # small tolerance to account for floating-point errors
    q = abs.(k)
    q .= ifelse.(abs.(q .- 2π) .< tolerance, 0.0, ifelse.(q .> π, -(q .- 2π), q))
    q = abs.(q)
    q .+= 1e-4
    q .= ifelse.(q .> π, q .- 1e-4, q)
    return q
end

function get_kvec(i, j, k)
    temp = BZ * [i / nx - 0.5, j / ny - 0.5, k / nz - 0.5]
    return temp[1:dim]
end

function get_kpts()
    kpts_list = Vector{Vector{Float64}}()
    for i in 1:nx, j in 1:ny, k in 1:nz
        kvec = get_kvec(i, j, k)
        e = epsilon(1, kvec)
        if abs(e - mu) > wD
            continue
        end
        push!(kpts_list, kvec)
    end
    return kpts_list
end

function get_bandwidth()
    maxval = -1000
    minval = 1000
    for i in 1:200, j in 1:200, k in 1:200
        kvec = get_kvec(i, j, k)
        eps = epsilon(1, kvec)
        maxval = max(maxval, eps)
        minval = min(minval, eps)
    end
    return maxval - minval
end

function paper_V(w)
    lambda = 2.0
    return lambda * 0.01 / (imag(w)^2 + 0.01)
end

function paper2_V(w)
    n = beta * imag(w) / (2 * pi)
    lambda = 2.0
    v = beta * wD / (2 * pi)
    return -lambda * v^2 / (n^2 + v^2)
end

function BCS_V(k, w)
    e = epsilon(1, k) - mu
    #e = 0.0
    lambda = -5.0
    if abs(e) < wD
        return lambda
    end
    return 0.0
end

function fill_V_arr(iw, frequency_dependence, V_obj)
    Vw_arr = Array{ComplexF32}(undef, length(iw), nx, ny, nz)
    for i in 1:nx, j in 1:ny, k in 1:nz, l in 1:length(iw)
        q = get_kvec(i, j, k)
        w = iw[l]
        Vw_arr[l, i, j, k] = V_obj(q, imag(w))
    end
    #Vw_arr .= Vw_arr .+ paper2_V.(iw)
    #Vw_arr .= -5.0
    if frequency_dependence
        return kw_to_rtau(Vw_arr, 'B', mesh)
    else
        return fft(Vw_arr)
    end
end

function initialize_phi_Z_chi!(phi_arr, Z_arr, chi_arr, iw)
    @threads for i in 1:nw
        @inbounds for j in 1:nx
            @inbounds for k in 1:ny
                @inbounds for l in 1:nz
                    w = iw[i]
                    kvec = get_kvec(j, k, l)
                    #phi_arr[i, j, k, l] = 0.02 / (imag(w)^2/2.0 + 1)# * dwave
                    phi_arr[i, j, k, l] = 1.0
                    Z_arr[i, j, k, l] = 1.0
                    chi_arr[i, j, k, l] = 0.0
                end
            end
        end
    end
end

function condense_to_F_and_G!(phi, Z, chi, F_arr, G_arr, iw, sigma, mu, projection)
    @threads for l in 1:length(iw)
        @inbounds for i in 1:nx
            @inbounds for j in 1:ny
                @inbounds for k in 1:nz
                    w = iw[l]
                    e = epsilon(1, get_kvec(i, j, k))
                    #e = 0.0
                    phi_el = phi[l] * projection[i, j, k]
                    denom = get_denominator(imag(w), phi_el, Z[l, i, j, k], chi[l, i, j, k] + e - mu)
                    F_arr[l, i, j, k] = -phi_el / denom
                    G_arr[l, i, j, k] = -(w * Z[l, i, j, k] + e - mu + chi[l, i, j, k]) / denom
                end
            end
        end
    end
end

function condense_to_F_and_G!(phi, Z, chi, F_arr, G_arr, iw, sigma, mu, e_arr, frequency_dependence)
    F_arr .= -phi ./ (-(Z .* iw) .^ 2 .+ phi .^ 2 .+ (e_arr .+ chi .- mu) .^ 2)
    if frequency_dependence
        F_arr .= F_arr .* ifelse.(abs.(e_arr .- mu) .> wD, 0.0, 1.0)
    end
    G_arr .= -(iw .* Z .+ e_arr .+ chi .- mu) ./ (-(Z .* iw) .^ 2 .+ phi .^ 2 .+ (e_arr .+ chi .- mu) .^ 2)
end

function update!(F, G, V_arr, phi, Z, chi, iw, sigma, projection)
    F_rt = kw_to_rtau(F, 'F', mesh)
    G_rt = kw_to_rtau(G, 'F', mesh)

    phit = V_arr .* F_rt / 2
    sigmat = -V_arr .* G_rt / 2

    phi_full = rtau_to_kw(phit, 'F', mesh)
    sigma .= rtau_to_kw(sigmat, 'F', mesh)

    phi .= phi_full[:, 1, 1, 1] / projection[1, 1, 1]
    fill_Z_chi!(iw, sigma, Z, chi)
end

function update!(F, G, V_arr, phi, Z, chi, iw, sigma)
    F_rt = kw_to_rtau(F, 'F', mesh)
    G_rt = kw_to_rtau(G, 'F', mesh)

    phit = V_arr .* F_rt / 2
    sigmat = -V_arr .* G_rt / 2

    phi .= rtau_to_kw(phit, 'F', mesh)
    sigma .= rtau_to_kw(sigmat, 'F', mesh)

    fill_Z_chi!(iw, sigma, Z, chi)
end

function update_phi!(F, V_arr, phi)
    F_rt = kw_to_rtau(F, 'F', mesh)
    phit = V_arr .* F_rt
    phi .= rtau_to_kw(phit, 'F', mesh)
end

function update_sigma!(G, V_arr, sigma)
    G_rt = kw_to_rtau(G, 'F', mesh)
    sigmat = -V_arr .* G_rt
    sigma .= rtau_to_kw(sigmat, 'F', mesh)
end

function update2!(F, G, V_arr, phi, Z, chi, iw, sigma)
    update_phi!(F, V_arr, phi)
    update_sigma!(G, V_arr, sigma)
    fill_Z_chi!(iw, sigma, Z, chi)
end

function update_fourier_transform!(F, G, V_arr, phi, Z, chi, iw, sigma, projection)
    F_rt = fft(F)
    G_rt = fft(G)

    phit = V_arr .* F_rt
    sigmat = -V_arr .* G_rt

    phi_full = fftshift(ifft(phit)) / (beta * nx * ny * nz)
    sigma .= fftshift(ifft(sigmat)) / (beta * nx * ny * nz)


    ind = argmax(abs.(phi_full))
    inds = Tuple(CartesianIndices(phi_full)[ind])
    phi .= phi_full[:, inds[2], inds[3], inds[4]] / projection[inds[2], inds[3], inds[4]]
    fill_Z_chi!(iw, sigma, Z, chi)
end

function update_fourier_transform!(F, G, V_arr, phi, Z, chi, iw, sigma)
    F_rt = fft(F)
    G_rt = fft(G)

    phit = V_arr .* F_rt
    sigmat = -V_arr .* G_rt

    phi .= fftshift(ifft(phit)) / (beta * nx * ny * nz)
    sigma .= fftshift(ifft(sigmat)) / (beta * nx * ny * nz)


    #ind = argmax(abs.(phi_full))
    #inds = Tuple(CartesianIndices(phi_full)[ind])
    #phi .= phi_full[:, inds[2], inds[3], inds[4]] / projection[inds[2], inds[3], inds[4]]
    fill_Z_chi!(iw, sigma, Z, chi)
end

function fill_Z_chi!(iw, sigma, Z, chi)
    #Z .= 1.0
    #chi .= 0.0
    #return
    Z .= 1.0 .- 0.5 .* (sigma .- sigma[end:-1:1, :, :, :]) ./ iw
    chi .= 0.5 .* (sigma .+ sigma[end:-1:1, :, :, :])
end


function find_chemical_potential(mu_initial, phi_arr, Z_arr, chi_arr, F_arr, G_arr, iw, sigma, target_Ne, projection)
    return mu_initial
    f(mu) = begin
        condense_to_F_and_G!(phi_arr, Z_arr, chi_arr, F_arr, G_arr, iw, sigma, mu, projection)
        get_number_of_electrons(G_arr) - target_Ne
    end
    return fzero(f, mu_initial; atol=1e-3, rtol=1e-3)
end

function get_number_of_electrons(G)
    Ne = 0.0
    for i in 1:nx, j in 1:ny, k in 1:nz
        Gsum = sum(real.(G[:, i, j, k]))
        Ne += 1 + 2 * Gsum / beta
    end
    return round(Ne / nk, digits=6)
end

function eliashberg_global()
    println("Beginning Eliashberg")
    IR_tol = 1e-5
    scf_tol = 1e-4

    #vertex = Firefly.Field_C(outdir * prefix * "_vertex.dat")
    vertex = Firefly.Vertex()

    frequency_dependence = true
    if frequency_dependence
        println("Creating IRMesh")
        global mesh = IR_Mesh(IR_tol)
        global fnw = mesh.fnw
        global fntau = mesh.fntau
        global bnw = mesh.bnw
        global bntau = mesh.bntau
        println("IRMesh created")
        iw = Array{ComplexF32}(undef, fnw)
        iv = Array{ComplexF32}(undef, bnw)
        for i in 1:fnw
            iw[i] = valueim(mesh.IR_basis_set.smpl_wn_f.sampling_points[i], beta)
        end
        for i in 1:bnw
            iv[i] = valueim(mesh.IR_basis_set.smpl_wn_b.sampling_points[i], beta)
        end

        println("Filled iw and iv")
        println("fnw, bnw: ", fnw, " ", bnw)
        println("Min & Max iw: ", minimum(imag.(iw)), " ", maximum(imag.(iw)))
        global nw = length(iw)
    else
        iw = Array{ComplexF32}(undef, nw)
        iv = Array{ComplexF32}(undef, nw)
        for i in 1:nw
            n = -nw / 2 + i - 1
            iw[i] = Complex(0.0, (2 * n - 1) * pi / beta)
            iv[i] = Complex(0.0, 2 * n * pi / beta)
        end
    end

    println("Initializing phi, Z, and chi")
    phi_arr = ones(Complex{Float32}, nw, nx, ny, nz) #* 1e-3
    #phi_arr .= 0.001 ./ (imag.(iw).^2 + 1)
    Z_arr = ones(Complex{Float32}, nw, nx, ny, nz)
    chi_arr = zeros(Complex{Float32}, nw, nx, ny, nz)

    println("Initializing epsilon")
    e_arr = zeros(Float32, nx, ny, nz)
    for i in 1:nx, j in 1:ny, k in 1:nz
        e_arr[i, j, k] = epsilon(1, get_kvec(i, j, k))
    end
    e_4d = reshape(e_arr, 1, nx, ny, nz)
    sigma = zeros(Complex{Float32}, nw, nx, ny, nz)
    println("Initializing V")
    V_arr = fill_V_arr(iv, frequency_dependence, vertex)

    println("Initializing F and G")
    F_arr = Array{ComplexF32}(undef, nw, nx, ny, nz)
    G_arr = Array{ComplexF32}(undef, nw, nx, ny, nz)
    #condense_to_F_and_G!(phi_arr, Z_arr, chi_arr, F_arr, G_arr, iw, sigma, mu, projection)
    condense_to_F_and_G!(phi_arr, Z_arr, chi_arr, F_arr, G_arr, iw, sigma, mu, e_4d, frequency_dependence)

    #Ne = round(get_number_of_electrons(G_arr), digits=6)

    phierr = 0.0
    prev_phi_arr = copy(phi_arr)
    println("Starting Convergence Loop")
    iterations = 1000
    #iterations = 1
    for i in 1:iterations
        if frequency_dependence
            #update!(F_arr, G_arr, V_arr, phi_arr, Z_arr, chi_arr, iw, sigma, projection)
            update2!(F_arr, G_arr, V_arr, phi_arr, Z_arr, chi_arr, iw, sigma)
        else
            #update_fourier_transform!(F_arr, G_arr, V_arr, phi_arr, Z_arr, chi_arr, iw, sigma, projection)
            update_fourier_transform!(F_arr, G_arr, V_arr, phi_arr, Z_arr, chi_arr, iw, sigma)
        end

        phierr = minimum((maximum(abs.(real.(phi_arr - prev_phi_arr))), maximum(abs.(real.(phi_arr + prev_phi_arr)))))
        max_phi = round(maximum(abs.(phi_arr)), digits=6)
        print("Iteration $i: MaxPhi = $max_phi              Error = $phierr           \r")

        if abs(phierr) < scf_tol || max_phi < 1e-10
            println("Converged after ", i, " iterations                                                                   ")
            break
        end
        #p = plot(real.(phi_arr[:, 1, 1, 1]))
        #display(p)
        #readline()
        prev_phi_arr = copy(phi_arr)
        #global mu = find_chemical_potential(mu, phi_arr, Z_arr, chi_arr, F_arr, G_arr, iw, sigma, Ne, projection)
        #condense_to_F_and_G!(phi_arr, Z_arr, chi_arr, F_arr, G_arr, iw, sigma, mu, projection)
        condense_to_F_and_G!(phi_arr, Z_arr, chi_arr, F_arr, G_arr, iw, sigma, mu, e_4d, frequency_dependence)
    end

    max_phi = maximum(abs.(phi_arr))
    max_Z = maximum(abs.(Z_arr))

    println("Error: ", phierr, "            ")
    println("Max phi: ", max_phi)
    println("Max Z: ", max_Z)
    return
    reordered_phi = Array{ComplexF32}(undef, nx, ny, nz, nw)
    reordered_Z = Array{ComplexF32}(undef, nx, ny, nz, nw)
    points = Matrix{Float32}(undef, nx * ny * nz * nw, dim + 1)
    @threads for i in 1:nw
        for j in 1:nx, k in 1:ny, l in 1:nz
            reordered_phi[j, k, l, i] = phi_arr[i, j, k, l]
            reordered_Z[j, k, l, i] = Z_arr[i, j, k, l]
            w = iw[i]
            kvec = get_kvec(j, k, l)
            idx = (i - 1) * nx * ny * nz + (j - 1) * ny * nz + (k - 1) * nz + l
            points[idx, 1:dim] .= kvec
            points[idx, dim+1] = imag(w)
        end
    end
    println("Reordered phi and Z")

    npts = Array{Integer}(undef, nw)
    for i in 1:nw
        w = iw[i]
        n = trunc((beta * imag(w) / pi - 1) / 2)
        npts[i] = n
    end
    println("Calculated npts")

    Z_index = argmax(abs.(Z_arr))
    Z_indices = Tuple(CartesianIndices(Z_arr)[Z_index])
    Z_w = Z_arr[:, Z_indices[2], Z_indices[3], Z_indices[4]]

    phi_average = average_over_k(real.(reordered_phi))
    Z_average = average_over_k(real.(reordered_Z))
    phi_average_imag = average_over_k(imag.(reordered_phi))
    println("Calculated averages")

    p = plot(npts, real.(phi_average), label="Phi", xlims=(-500, 500))
    #plot!(p, phi_average_imag, label="Imag Phi")
    display(p)
    readline()
    p = plot(npts, real.(Z_average), label="Z", xlims=(-500, 500))
    display(p)
    readline()

    println("Saving phi and Z")
    #save_data(prefix*"_phi.dat", points, reordered_phi, dim, true, true, false)
    #save_data(prefix*"_Z.dat", points, reordered_Z, dim, true, true, false)
end

function get_denominator(w1, phi_el, Z_el, eps)
    return abs2(w1 * real(Z_el)) + abs2(eps) + abs2(real(phi_el))
end

function average_over_k(arr)
    #return Float32.(arr[1, 1, 1, :])
    average = zeros(Float32, length(arr[1, 1, 1, :]))
    for i in nx, j in ny, k in nz
        average .+= (arr[i, j, k, :])
    end
    return average
end

function fill_gpu!(F_arr, V_arr, e_arr, phi, pts, nx, ny, beta, BZ, mu, w_arr, F_func, V_func, e_func)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x  # Compute global index
    if i ≤ pts
        for j in 1:nx, k in 1:ny
            n = -pts / 2 + i - 1
            iwn = Complex(0.0, (2 * n + 1) * pi / beta)
            iwm = Complex(0.0, (2 * n) * pi / beta)
            w_arr[i] = iwn
            kx = BZ[1, 1] * (j / nx - 0.5) + BZ[1, 2] * (k / ny - 0.5) + BZ[1, 3]
            ky = BZ[2, 1] * (j / nx - 0.5) + BZ[2, 2] * (k / ny - 0.5) + BZ[2, 3]
            e_arr[j, k] = e_func(kx, ky, 0) - mu
            F_arr[i, j, k] = F_func(iwn, e_arr[j, k], phi[i])
            V_arr[i, j, k] = V_func(iwm, e_arr[j, k])
        end
    end
end

function update_gpu!(F_arr, e_arr, phi, pts, nx, ny, w_arr, F_func)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x  # Compute global index
    if i ≤ pts
        for j in 1:nx, k in 1:ny
            F_arr[i, j, k] = F_func(w_arr[i], e_arr[j, k], phi[i, 1, 1, 1])
        end
    end
end

function error_gpu!(phi, prev_phi, error, max, pts, nx, ny)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x  # Compute global index
    if i ≤ pts
        max_val = 0
        for j in 1:nx, k in 1:ny
            error[i] += abs(phi[i, j, k] - prev_phi[i, j, k])
            if phi[i, j, k] > max_val
                max_val = phi[i, j, k]
            end
            prev_phi[i, j, k] = phi[i, j, k]
        end
        max[i] = max_val
    end
end

function fftshift_gpu!(result::CuArray)
    for i in 1:size(result, 1)
        for j in 1:size(result, 2)
            result[i, j] = result[i, j] * (-1)^(i + j)
        end
    end
end

function VF_3sum_GPU(pts)
    println("Starting")
    if dim != 2
        println("Only 2D supported")
        return
    end
    nk = nx * ny

    lambda = 32.0
    wD = 0.2
    mu = 1.0

    F(iw, e, phi) = phi / (-iw^2 + phi^2 + e^2)
    V(iv, e) = abs(e) > wD ? 0.0 : lambda

    F_arr = CUDA.zeros(ComplexF32, pts, nx, ny)
    V_arr = CUDA.zeros(ComplexF32, pts, nx, ny)
    e_arr = CUDA.zeros(Float32, nx, ny)
    w_arr = CUDA.zeros(ComplexF32, pts)
    phi = CUDA.ones(Float32, pts) * 1e-3
    BZ_GPU = CuArray(BZ)
    println("Filling arrays")
    threads = 256
    blocks_per_grid = cld(pts, threads)
    CUDA.@cuda threads = threads blocks = blocks_per_grid fill_gpu!(F_arr, V_arr, e_arr, phi, pts, nx, ny, beta, BZ_GPU, mu, w_arr, F, V, epsilon_cuda)

    println("Starting Fourier Transforms")
    error = CUDA.ones(Float32, pts)
    max_arr = CUDA.zeros(Float32, pts)
    for i in 1:50
        F_rt = CUDA.CUFFT.fft(F_arr)
        V_rt = CUDA.CUFFT.fft(V_arr)
        phi_rt = CUDA.dot(F_rt, V_rt)
        ifft = CUDA.CUFFT.ifft(phi_rt)
        result = CUDA.fftshift(CUDA.CUFFT.ifft(phi_rt)) / (pts * nk)

        CUDA.@cuda threads = threads blocks = blocks_per_grid error_gpu(result, phi, error, max_arr, pts, nx, ny)
        error_total = sum(Array(error))
        max_total = max(Array(max_arr))
        println("Error total: ", error_total)
        println("Max val: ", max_total)

        if error_total < 1e-4
            break
        end
        CUDA.@cuda threads = threads blocks = blocks_per_grid update_gpu!(F_arr, e_arr, result, pts, nx, ny, w_arr, F)
    end

    phi_result = Array(result)
    plot_result = phi_result[:, 1, 1, 1]
    p = plot(plot_result)
    display(p)
    readline()
end

function VF_sum()
    pts = 10000
    basis = IR_Mesh(1e-10)
    fnw, bnw, fntau, bntau = basis.fnw, basis.bnw, basis.fntau, basis.bntau
    result = zeros(ComplexF32, pts)
    F_arr = zeros(ComplexF32, pts)
    V_arr = zeros(ComplexF32, pts)

    w_arr = Array{Float32}(undef, pts)
    iw_arr = Array{Float32}(undef, fnw)
    for i in 1:fnw
        iw_arr[i] = imag(valueim(basis.IR_basis_set.smpl_wn_f.sampling_points[i], beta))
    end

    F_ir = Array{ComplexF32}(undef, fnw)
    V_ir = Array{ComplexF32}(undef, bnw)
    phi = 1.0
    lambda = 32.0
    wD = 0.2
    e = 0.4

    F(iw) = phi / (-iw^2 + phi^2 + e^2)
    V(iv) = lambda * 2 * wD^2 / (wD^2 - iv^2)

    for i in 1:fnw
        iw = valueim(basis.IR_basis_set.smpl_wn_f.sampling_points[i], beta)
        F_ir[i] = F(iw)
    end
    for i in 1:bnw
        iv = valueim(basis.IR_basis_set.smpl_wn_b.sampling_points[i], beta)
        V_ir[i] = V(iv)
    end

    @threads for i in 1:pts
        for j in 1:pts
            n = -pts / 2 + i - 1
            m = -pts / 2 + j - 1
            iwn = Complex(0.0, (2 * n + 1) * pi / beta)
            w_arr[i] = imag(iwn)
            iwm = Complex(0.0, (2 * m + 1) * pi / beta)
            iv = Complex(0.0, 2 * n * pi / beta)
            result[i] += (V(iwn - iwm) * F(iwm)) / beta
            F_arr[i] = F(iwn)
            V_arr[i] = V(iv)
        end
    end

    F_rt = fft(F_arr)
    V_rt = fft(V_arr)
    phi_rt = F_rt .* V_rt
    result2 = fftshift(ifft(phi_rt) / pts)

    F_irt = wn_to_tau(basis, Fermionic(), F_ir)
    V_irt = wn_to_tau(basis, Bosonic(), V_ir)
    phi_irt = F_irt .* V_irt
    result3 = tau_to_wn(basis, Fermionic(), phi_irt)

    p = plot(w_arr, real.(result), label="Summed", xlims=(-20, 20))
    plot!(p, w_arr, real.(result2), label="Ft'd", xlims=(-20, 20))
    plot!(p, iw_arr, real.(result3), label="IR'd", xlims=(-20, 20))
    display(p)
    readline()
end

function VF_4sum()
    pts = 10000
    println("Starting")
    basis = IR_Mesh(1e-10)
    fnw, bnw, fntau, bntau = basis.fnw, basis.bnw, basis.fntau, basis.bntau
    nk = nx * ny * nz
    result = zeros(ComplexF32, pts)
    F_arr = zeros(ComplexF32, pts, nx, ny, nz)
    V_arr = zeros(ComplexF32, pts, nx, ny, nz)

    w_arr = Array{Float32}(undef, pts)
    iw_arr = Array{Float32}(undef, fnw)
    for i in 1:fnw
        iw_arr[i] = imag(valueim(basis.IR_basis_set.smpl_wn_f.sampling_points[i], beta))
    end

    F_ir = Array{ComplexF32}(undef, fnw, nx, ny, nz)
    V_ir = Array{ComplexF32}(undef, bnw, nx, ny, nz)
    phi = 1.0
    lambda = 32.0
    wD = 0.2

    F(iw, e) = phi / (-iw^2 + phi^2 + e^2)
    V(iv) = lambda * 2 * wD^2 / (wD^2 - iv^2)

    println("Filling IR arrays")
    for i in 1:fnw, j in 1:nx, k in 1:ny, l in 1:nz
        iw = valueim(basis.IR_basis_set.smpl_wn_f.sampling_points[i], beta)
        F_ir[i, j, k, l] = F(iw, epsilon(get_kvec(j, k, l)))
    end
    for i in 1:bnw, j in 1:nx, k in 1:ny, l in 1:nz
        iv = valueim(basis.IR_basis_set.smpl_wn_b.sampling_points[i], beta)
        V_ir[i, j, k, l] = V(iv)
    end

    println("Starting sum")
    for i in 1:pts
        print("Percent complete: ", round(i / pts * 100, digits=2), "%                         \r")
        @threads for j in 1:pts
            n = -pts / 2 + i - 1
            m = -pts / 2 + j - 1
            iwn = Complex(0.0, (2 * n + 1) * pi / beta)
            w_arr[i] = imag(iwn)
            iwm = Complex(0.0, (2 * m + 1) * pi / beta)
            iv = Complex(0.0, 2 * n * pi / beta)
            @inbounds for a in 1:nx
                @inbounds for b in 1:ny
                    @inbounds for c in 1:nz
                        em = epsilon(get_kvec(a, b, c))
                        result[i] += (V(iwn - iwm) * F(iwm, em)) / (beta * nk)
                        F_arr[i, a, b, c] = F(iwn, epsilon(get_kvec(a, b, c)))
                        V_arr[i, a, b, c] = V(iv)
                    end
                end
            end
        end
    end

    println("Starting Fourier Transforms")
    F_rt = fft(F_arr)
    V_rt = fft(V_arr)
    phi_rt = F_rt .* V_rt
    result2 = fftshift(ifft(phi_rt)) / (pts * nk)

    F_irt = kw_to_rtau(F_ir, 'F', basis)
    V_irt = kw_to_rtau(V_ir, 'B', basis)
    phi_irt = F_irt .* V_irt
    result3 = rtau_to_kw(phi_irt, 'F', basis)

    println("Printing max vals")
    println(maximum(real.(result)), " ", maximum(real.(result2)), " ", maximum(real.(result3)))


    println("Plotting")
    result_plot = result
    result2_plot = sum(result2, dims=(2, 3, 4))[:, 1, 1, 1] / (nx * ny * nz)
    result3_plot = sum(result3, dims=(2, 3, 4))[:, 1, 1, 1] / (nx * ny * nz)

    p = plot(w_arr, real.(result_plot), label="Summed", xlims=(-20, 20))
    plot!(p, w_arr, real.(result2_plot), label="Ft'd", xlims=(-20, 20))
    plot!(p, iw_arr, real.(result3_plot), label="IR'd", xlims=(-20, 20))
    display(p)
    readline()
end


@inline function F_mod(iw, e, phi)
    val = phi / (-iw^2 + phi^2 + e^2)
    if abs(e) < 0.01
        return val
    else
        return val * 0.0
    end
    #return phi / (-iw^2 + phi^2 + e^2)
end
function fill_F_arr_gpu!(F_arr, V_arr, phi, pts, nx, ny, beta, BZ, mu, F_func, V_func, e_func, wD, lambda)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x  # Compute global index
    if i ≤ pts
        for j in 1:nx, k in 1:ny
            n = -pts / 2 + i - 1
            iwn = Complex(0.0, (2 * n + 1) * pi / beta)
            iwm = Complex(0.0, (2 * n) * pi / beta)
            kx = BZ[1, 1] * (j / nx - 0.5) + BZ[1, 2] * (k / ny - 0.5) + BZ[1, 3]
            ky = BZ[2, 1] * (j / nx - 0.5) + BZ[2, 2] * (k / ny - 0.5) + BZ[2, 3]
            e = e_func(kx, ky, 0) - mu
            F_arr[i, j, k] = F_func(iwn, e, phi[i])
            V_arr[i, j, k] = V_func(iwm, e, wD, lambda)
        end
    end
end

function get_transformed_V()
    V = zeros(Float32, nw, nx, ny)
    #V .= V .+ ifelse.(abs.(e_3d) .> wD, 0.0, 10.0)
    #V .= V .+ V_gpu.(iv, e_3d, 0.2, 10.0)
    V .= V .+ 1.0
    return fft(V)
end

function energy_fastest()
    nk = nx * ny * nz
    println("nk: ", nk)
    phi = ones(Float64, nx, ny, nz)
    println("Filling epsilon array")
    e = zeros(Float64, nx, ny, nz)
    kpts = Array{Vector{Float64}}(undef, nx, ny, nz)
    for i in 1:nx, j in 1:ny, k in 1:nz
        kpts[i, j, k] = get_kvec(i, j, k)
    end
    e = epsilon.(1, kpts) .- mu
    println("Filling V array")
    V = ones(Float64, nx, ny, nz) .* 1.0 #.* ifelse.(abs.(e) .> 0.1, 0.0, 1.0)
    V_rt = fft(V)
    println("Starting iterations")
    dos = Firefly.Field_R(outdir * prefix * "_DOS.dat")

    initial = dos(mu) * asinh(wD)
    println("Expected initial: ", initial)
    #sumval = 0.00018
    iterations = 1
    diff = 0
    for i in 1:iterations
        F = phi ./ (2 .* (phi .^ 2 .+ e .^ 2) .^ 0.5) .* ifelse.(abs.(e) .> wD, 0.0, 1.0)
        F_rt = fft(F)
        phi_rt = F_rt .* V_rt
        result = fftshift(ifft(phi_rt)) / (nk)
        phi = result
        diff = abs(diff - maximum(abs.(phi)))
        if (i - 1) % 5 == 0
            println(i, ": Max val: ", maximum(abs.(phi)))
        end
    end
    println("Diff: ", diff)
    final = wD / sinh(1 / dos(mu))
    println("Expected final: ", final)
end

function func_gap_int(E, T)
    if abs(E) < 0.0001
        return 1 / (4 * T)
    end
    return tanh(E / (2 * T)) / (2 * E)
end


function energy_integral(D)
    dos = Firefly.Field_R(outdir * prefix * "_DOS.dat")
    emin = -wD
    emax = wD
    npts = 100000
    sum = 0
    T = cfg.Temperature
    for i in 1:npts
        e = emin + (emax - emin) * (i - 1) / (npts - 1)
        E = (e^2 + D^2)^0.5
        f = func_gap_int(E, T)
        sum += dos(e + mu) * f
        #println(func_gap_int(E, T))
        #sum += dos(e)
    end
    sum *= (emax - emin) / npts
    return sum
end

function energy_finite_T()
    D = 0.001
    dos = Firefly.Field_R(outdir * prefix * "_DOS.dat")
    sample_f = func_gap_int(D, cfg.Temperature)
    println("DOS = ", dos(mu))
    println("calculated = ", energy_integral(D))
    println("theory = ", dos(mu) * asinh(wD / D))
    println("debug_estimate = ", 2 * wD * dos(mu) * sample_f)
    #for i in 1:10
    #    D = i * D_init
    #    println("D=", D, ", 1=", energy_integral(D))
    #end
end

function fastest()
    nk = nx * ny * nz
    println("nk: ", nk)
    phi = ones(Float32, nw)
    n = [-nw / 2 + i - 1 for i in 1:nw]
    iw = [Complex(0.0, (2 * n + 1) * pi / beta) for n in n]
    iv = [Complex(0.0, 2 * n * pi / beta) for n in n]
    println("Filling epsilon array")
    e = zeros(Float32, nx, ny)
    for i in 1:nx, j in 1:ny
        kvec = get_kvec(i, j, 100)
        e[i, j] = epsilon(1, kvec) - mu
    end
    println("Filling V array")
    iw_3d = reshape(iw, (nw, 1, 1))
    e_3d = reshape(e, (1, nx, ny))
    phi_3d = reshape(phi, (nw, 1, 1))
    V_rt = get_transformed_V()
    println("Starting iterations")

    dos = Firefly.Field_R(outdir * prefix * "_DOS.dat")
    initial = dos(mu) * asinh(wD)
    println("Expected initial: ", initial)
    iterations = 8
    for i in 1:iterations
        F = F_mod.(iw_3d, e_3d, phi_3d)
        F_rt = fft(F)
        phi_rt = F_rt .* V_rt
        result = fftshift(ifft(phi_rt)) / (beta * nk)
        println("Max val: ", maximum(abs.(result)))
        phi_3d = result
    end
    final = wD / sinh(1 / dos(mu))
    println("Expected final: ", final)
end

@inline function V_gpu(iv, e, wD, lambda)
    return abs(e) > wD ? 0.0 : lambda
    #return lambda
    #return lambda * 2 * wD^2 / (wD^2 - iv^2)
end

@inline function F_gpu(iw, e, phi)
    return phi / (-iw^2 + phi^2 + e^2)
end

@inline function epsilon_cuda(kx, ky, kz)
    return -2 * (cos(kx) + cos(ky))
    #return kx^2 + ky^2
end


function compute_gpu_sum!(result, pts, nx, ny, nz, beta, phi, wD, lambda, w_arr_CUDA, BZ, mu)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x

    if i ≤ pts
        iwn = ComplexF32(w_arr_CUDA[i])
        sum = 0

        for j in 1:pts, a in 1:nx, b in 1:ny, c in 1:nz
            iwm = ComplexF32(w_arr_CUDA[j])
            dw = iwn - iwm
            kx = BZ[1, 1] * (a / nx - 0.5) + BZ[1, 2] * (b / ny - 0.5) + BZ[1, 3] * (c / nz - 0.5)
            ky = BZ[2, 1] * (a / nx - 0.5) + BZ[2, 2] * (b / ny - 0.5) + BZ[2, 3] * (c / nz - 0.5)
            kz = BZ[3, 1] * (a / nx - 0.5) + BZ[3, 2] * (b / ny - 0.5) + BZ[3, 3] * (c / nz - 0.5)
            em = epsilon_cuda(kx, ky, kz) - mu
            sum += (V_gpu(dw, em, wD, lambda) * F_gpu(iwm, em, phi[j])) / (beta * nx * ny * nz)
        end
        result[i] = sum
    end
    return nothing
end

function energy_2sum()
    println("Starting")
    nk = nx * ny * nz
    println("nk: ", nk)

    lambda = 1.0
    wD = cfg.cutoff_energy
    dos = Firefly.Field_R(outdir * prefix * "_DOS.dat")

    println("Starting CPU sum")
    result::Float64 = 0.0
    for i in 1:nx, j in 1:ny, k in 1:nz
        e = epsilon(1, get_kvec(i, j, k)) - mu
        if abs(e) < wD
            result += lambda
        end
    end
    println("Result: ", result / nk)
    println("Expected: ", 2 * wD * dos(mu))
    return nothing
end

function VF_4sum_GPU(phi, pts, phi_ir)
    println("Starting")
    nk = nx * ny * nz
    println("nk: ", nk)

    lambda = 32.0
    wD = 0.2
    mu = 1.0

    F(iw, e, phi) = phi / (-iw^2 + phi^2 + e^2)
    V(iv, e) = abs(e) > wD ? 0.0 : lambda

    result = CUDA.zeros(ComplexF32, pts)
    F_arr = CUDA.ones(ComplexF32, pts, nx, ny)
    V_arr = CUDA.ones(ComplexF32, pts, nx, ny)
    w_arr = zeros(ComplexF32, pts)
    w_array = Array{ComplexF32}(undef, pts)
    println("Filling arrays")
    @threads for i in 1:pts
        n = -pts / 2 + i - 1
        w_arr[i] = Complex(0.0, (2 * n - 1) * pi / beta)
        w_array[i] = Complex(0.0, (2 * n - 1) * pi / beta)
    end
    w_arr_CUDA = CuArray(w_array)
    BZ_CUDA = CuArray(BZ)

    phi_CUDA = CuArray(phi)
    println("Starting GPU sum")

    CUDA.@sync @cuda threads = 512 blocks = cld(pts, 512) compute_gpu_sum!(result, pts, nx, ny, nz, beta, phi_CUDA, wD, lambda, w_arr_CUDA, BZ_CUDA, mu)
    cpu_result = Array(result)
    println("Max result: ", maximum(abs.(cpu_result)))
    return cpu_result, ones(ComplexF32, 3)
end

function gpu_eliashberg()
    pts = nw
    phi = ones(ComplexF32, pts)# * 1e-3
    #basis = FiniteTempBasisSet(beta, Float32(1.0), IR_tol)
    #fnw, bnw, fntau, bntau = basis.fnw, basis.bnw, basis.fntau, basis.bntau
    phi_ir = ones(ComplexF32, 3)
    prev_phi = ones(ComplexF32, pts)
    iterations = 1
    for i in 1:iterations
        prev_phi = copy(phi)
        println("Iteration: ", i)
        phi, phi_ir = Eliashberg.VF_4sum_GPU(phi, pts, phi_ir)
        phi_error = maximum(abs.(phi - prev_phi))
        if phi_error < 1e-6
            println("Converged")
            break
        end
    end
    p = plot(real.(phi))
    display(p)
    readline()
end

function element_convolution_F(A, B, mesh)
    A_t = wn_to_tau(mesh, Bosonic(), A)
    B_t = wn_to_tau(mesh, Fermionic(), B)
    temp = A_t .* B_t
    return tau_to_wn(mesh, Fermionic(), temp)
end

function element_convolution_B(A, B, mesh)
    A_t = wn_to_tau(mesh, Fermionic(), A)
    B_t = wn_to_tau(mesh, Fermionic(), B)
    temp = A_t .* B_t
    return tau_to_wn(mesh, Bosonic(), temp)
end

function element_convolution_const(A, B)
    A_t = fft(A)
    B_t = fft(B)
    temp = A_t .* reverse(B_t)
    return fftshift(ifft(temp)) / beta
end

function get_iw_iv_mesh(particle_type)
    method = 1
    if particle_type == "F" || particle_type == "Fermionic"
        method = 2
    elseif particle_type == "B" || particle_type == "Bosonic"
        method = 3
    end
    fnw = nw
    bnw = nw
    if method > 1
        mesh = IR_Mesh()
        fnw = mesh.fnw
        bnw = mesh.bnw
        iw = Array{ComplexF64}(undef, fnw)
        iv = Array{ComplexF64}(undef, bnw)
        for i in 1:fnw
            iw[i] = valueim(mesh.IR_basis_set.smpl_wn_f.sampling_points[i], beta)
        end
        for i in 1:bnw
            iv[i] = valueim(mesh.IR_basis_set.smpl_wn_b.sampling_points[i], beta)
        end
    else
        mesh = 1
        iw = Array{ComplexF64}(undef, nw)
        iv = Array{ComplexF64}(undef, nw)
        for i in 1:nw
            n = -nw / 2 + i - 1
            iw[i] = Complex(0.0, (2 * n - 1) * pi / beta)
            iv[i] = Complex(0.0, 2 * n * pi / beta)
        end
    end
    return iw, iv, mesh
end


# particle_type is of resultant particle. So if a boson is returned, it is bosonic
function conv_sum(particle_type, kvecs, func1, func2)
    iw, iv, mesh = get_iw_iv_mesh(particle_type)
    nk = length(kvecs)
    fnw = length(iw)
    bnw = length(iv)
    println(fnw)
    println(bnw)

    println("2)")
    P = [zeros(ComplexF64, bnw) for _ in 1:nk]
    Gw = zeros(ComplexF64, fnw)
    Fw = zeros(ComplexF64, fnw)
    if particle_type == "F" || particle_type == "Fermionic"
        P = [zeros(ComplexF64, fnw) for _ in 1:nk]
        #Gw = zeros(ComplexF64, fnw)
        Fw = zeros(ComplexF64, bnw)
    end

    println("3)")
    for j in 1:nk
        k1 = kvecs[j]
        #Gw .= 1 ./ (iw .- epsilon(1, k1) .+ mu)
        if particle_type == "F" || particle_type == "Fermionic"
            Gw .= func1(k1, iw)
        end

        println("4")
        sum = zeros(ComplexF64, bnw)
        for k in 1:nk
            q = kvecs[k] - k1
            #Fw .= 1 ./ (iw .- epsilon(1, q) .+ mu)
            if particle_type == "B" || particle_type == "Bosonic"
                Fw .= func2(q, iv)
            end

            println("5")
            if particle_type == "B" || particle_type == "Bosonic"
                sum += element_convolution_B(Fw, Gw, mesh)
            elseif particle_type == "F" || particle_type == "Fermionic"
                sum += element_convolution_F(Fw, Gw, mesh)
            else
                sum += element_convolution_const(Fw, Gw)
            end
        end
        print(j, " out of ", nk, "\r")
        P[j] = sum
    end
    return P
end

function conv_sum_eli(kvecs, interaction, phi, Z, chi, iw, iv, mesh)
    nk = length(kvecs)
    fnw = length(iw)
    bnw = length(iv)

    Gw = zeros(ComplexF64, fnw)
    P = [zeros(ComplexF64, fnw) for _ in 1:nk]
    Fw = zeros(ComplexF64, bnw)

    for j in 1:nk
        k1 = kvecs[j]
        Gw .= F_eliashberg(phi[j], Z[j], chi[j], iw, epsilon(1, k1))

        sum = zeros(ComplexF64, fnw)
        for k in 1:nk
            q = kvecs[k] - k1
            Fw .= get_V(interaction, q, iw, iv)

            sum += element_convolution_F(Fw, Gw, mesh)
        end
        P[j] = sum
    end
    return P
end

function conv_sum_const(F::Array{Vector{ComplexF64}}, G::Vector{Vector{ComplexF64}})
    nk = length(G)
    P = [zeros(ComplexF64, nw) for _ in 1:nk]  # P is Vector{Vector{ComplexF64}}
    for j in 1:nk
        Gw = G[j]
        sum = zeros(ComplexF64, nw)
        for k in 1:nk
            Fw = F[j, k]
            sum += element_convolution_const(Fw, Gw)
        end
        P[j] = sum
    end
    return P
end

function partial_convolution(V, phi, iw, iv, kvecs)
    nk = length(phi)
    fnw = length(iw)
    bnw = length(iv)
    P = [zeros(ComplexF64, fnw) for _ in 1:nk]  # P is Vector{Vector{ComplexF64}}
    Vw = zeros(ComplexF64, bnw)
    Fw = zeros(ComplexF64, fnw)

    println("Number of Julia threads: ", Threads.nthreads())
    for i in 1:nk
        phiw = phi[i]
        e = epsilon(1, kvecs[i]) - mu
        Fw .= abs(e) < wD ? phiw ./ (-iw .^ 2 .+ phiw .^ 2 .+ e^2) : 0
        k1 = kvecs[i]
        @threads for j in 1:nk
            dk = kvecs[j] - k1
            Vw .= ones(ComplexF64, bnw)
            #Vw .= [V(dk, iv[k]) for k in 1:bnw]
            X = element_convolution_const(Vw, Fw)
            P[j] += X
        end
    end
    return P / (nx * ny * nz)
end

function bcs_convsum()
    iw = Array{ComplexF32}(undef, nw)
    iv = Array{ComplexF32}(undef, nw)
    for i in 1:nw
        n = -nw / 2 + i - 1
        iw[i] = Complex(0.0, (2 * n - 1) * pi / beta)
        iv[i] = Complex(0.0, 2 * n * pi / beta)
    end

    kvecs = Vector{Vector{Float64}}(undef, nx * ny * nz)
    klen = length(kvecs)
    for i in 1:nx, j in 1:ny
        kvecs[(i-1)*ny+j] = get_kvec(i, j, 1)
    end
    phi = [ones(ComplexF64, nw) for _ in 1:klen]
    #V = Firefly.Vertex()
    V = ones(ComplexF64, nw)

    P = partial_convolution(V, phi, iw, iv, kvecs)
    println(typeof(P))
    println(maximum([maximum(abs.(real.(p))) for p in P]))
    println(P[0][0], P[14][15])
end

function F_eliashberg(phi, Z, chi, iw, e)
    theta = phi .^ 2 .+ (Z .* iw) .^ 2 + (chi .+ e) .^ 2
    return phi ./ theta
end

function get_V(V, k::Vector{Float64}, w, v)
    l = length(imag.(v))
    temp = Vector{ComplexF64}(undef, l)
    for i in 1:l
        val = 0
        if abs(w[i]) < wD && abs(v[i]) < wD
            val = 1
        end
        temp[i] = val
    end
    return temp
    return [V(k, imag(w[i])) for i in 1:l]
end

function get_max_phi(phi)
    max_real = -Inf
    for vec in phi
        for z in vec
            r = real(z)
            if r > max_real
                max_real = r
            end
        end
    end
    return max_real
end

function get_min_phi(phi)
    min_real = Inf
    for vec in phi
        for z in vec
            r = real(z)
            if r < min_real
                min_real = r
            end
        end
    end
    return min_real
end

function eliashberg_convsum()
    interaction = Firefly.Vertex()
    kvecs = Vector{Vector{Float64}}(undef, nx * ny)
    for i in 1:nx, j in 1:ny
        kvecs[(i-1)*ny+j] = get_kvec(i, j, 1)
    end
    nk = length(kvecs)

    iw, iv, mesh = get_iw_iv_mesh("F")
    fnw = length(iw)
    phi = [ones(ComplexF64, fnw) for _ in 1:nk]
    Z = [ones(ComplexF64, fnw) for _ in 1:nk]
    chi = [zeros(ComplexF64, fnw) for _ in 1:nk]

    iterations = 1
    for _ in 1:iterations
        phi = conv_sum_eli(kvecs, interaction, phi, Z, chi, iw, iv, mesh) / (nk)
        println("Max_Phi: ", get_max_phi(phi))
    end
    println("Max phi: ", get_max_phi(phi))
    println("Min phi: ", get_min_phi(phi))
end

function greens_function(k, iw)
    e = epsilon(1, k) - mu
    return 1 ./ (iw .- e)
end

function chi_convsum()
    kvecs = Vector{Vector{Float64}}(undef, nx * ny)
    for i in 1:nx, j in 1:ny
        kvecs[(i-1)*ny+j] = get_kvec(i, j, 1)
    end

    iw, iv, mesh = get_iw_iv_mesh("")
    P = conv_sum("", kvecs, greens_function, greens_function) / (nx * ny * nz)
    println("shape of P: ", size(P), " ", length(P), " ", length(P[1]))
    println(length(kvecs), " ", length(iv))
    open("convsum.dat", "w") do io
        println(io, "kx             ky             w             Re(f)            Im(f)")
        for i in 1:length(P), j in 1:length(P[1])
            kvec = kvecs[i]
            w = iv[j]
            temp = P[i][j]
            println(io,
                @sprintf("%.6f", kvec[1]), "    ",
                @sprintf("%.6f", kvec[2]), "    ",
                @sprintf("%.6f", imag(w)), "    ",
                @sprintf("%.6f", real(temp)), "    ",
                @sprintf("%.6f", imag(temp))
            )
        end
    end
end

function chi_conv_test()
    println("1")
    kvecs = Vector{Vector{Float64}}(undef, nx * ny)
    for i in 1:nx, j in 1:ny
        kvecs[(i-1)*ny+j] = get_kvec(i, j, 1)
    end
    println("2")
    iw, iv, mesh = get_iw_iv_mesh("")
    println("3")
    println(size(iw))
    fnw = length(iw)
    gkio = Array{ComplexF64}(undef, fnw, nx, ny, nz)
    for ix in 1:nx, iy in 1:ny, iw in 1:basis.fnw, iz in 1:nz
        iv::ComplexF64 = valueim(basis.IR_basis_set.smpl_wn_f.sampling_points[iw], beta)
        e = epsilon(1, kvec[(i-1)*ny+j])
        gkio[iw, ix, iy, iz] = 1.0 / (iv - e + mu)
    end
    println("4")
    grit = fft(gkio)
    println("5")
    crit = Array{ComplexF64}(undef, nw, nx, ny, nz)
    for ix in 1:nx, iy in 1:ny, iz in 1:nz, it in 1:basis.bntau
        crit[it, ix, iy, iz] = grit[it, ix, iy, iz] * grit[basis.bntau-it+1, ix, iy, iz]
    end
    println("6")
    ckio = fftshift(ifft(crit)) / (beta * nx * ny * nz)
    println("7")
    open("chi.dat", "w") do io
        println(io, "kx             ky             w             Re(f)            Im(f)")
        for i in 1:nw, j in 1:nx, k in 1:ny, l in 1:nz
            kvec = kvecs[(i-1)*ny+j]
            w = iv[i]
            temp = ckio[i, j, k, l]
            println(io,
                @sprintf("%.6f", kvec[1]), "    ",
                @sprintf("%.6f", kvec[2]), "    ",
                @sprintf("%.6f", imag(w)), "    ",
                @sprintf("%.6f", real(temp)), "    ",
                @sprintf("%.6f", imag(temp))
            )
        end
    end
end

function square_well_test()
    iw, iv, mesh = get_iw_iv_mesh("F")
    bnw, fnw = mesh.bnw, mesh.fnw
    V = Array{ComplexF64}(undef, bnw)
    F = Array{ComplexF64}(undef, fnw)

    a = 0.5
    w = 0.01
    for i in 1:bnw
        x = iv[i]
        V[i] = 0.5 * (tanh((x + a) / w) - tanh((x - a) / w))
    end
    for i in 1:fnw
        x = iw[i]
        F[i] = 0.5 * (tanh((x + a) / w) - tanh((x - a) / w))
    end

    Vt = wn_to_tau(mesh, Bosonic(), V)
    Ft = wn_to_tau(mesh, Fermionic(), F)
    temp = Vt .* Ft
    result = tau_to_wn(mesh, Fermionic(), temp)
    r_result = real.(result)
    println("max: ", maximum(r_result))
    println("min: ", minimum(r_result))
    println("avg: ", mean(r_result))

    V = Array{ComplexF64}(undef, nw)
    F = Array{ComplexF64}(undef, nw)

    a = 0.5
    w = 0.01
    for i in 1:nw
        n = -nw / 2 + i - 1
        iw = Complex(0.0, (2 * n - 1) * pi / beta)
        iv = Complex(0.0, (2 * n) * pi / beta)
        V[i] = 0.5 * (tanh((iv + a) / w) - tanh((iv - a) / w))
        F[i] = 0.5 * (tanh((iw + a) / w) - tanh((iw - a) / w))
    end


    Vt = fft(V)
    Ft = fft(F)
    result = fftshift(ifft(Vt .* Ft)) / beta
    r_result = real.(result)
    println("max: ", maximum(r_result))
    println("min: ", minimum(r_result))
    println("avg: ", mean(r_result))
end

function eliashberg_node()
    #energy_finite_T()
    #energy_2sum()
    #energy_fastest()
    #chi_conv_test()
    #chi_convsum()
    #square_well_test()
    #eliashberg_convsum()
    #fastest()
    eliashberg_global()
    println("k_mesh: ", cfg.k_mesh)
    return
end

end # module
