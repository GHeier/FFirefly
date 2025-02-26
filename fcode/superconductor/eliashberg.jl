module Eliashberg

t1 = time()
include("../objects/mesh.jl")
using .IRMesh

include("../objects/CondensedMatterField/julia_field.jl")
using .Field

using CUDA, FFTW
using Plots
using Roots
using SparseIR
import SparseIR: Statistics, value, valueim, MatsubaraSampling64F, TauSampling64
using Base.Threads, JLD
using LinearAlgebra, Printf, PyCall
np = pyimport("numpy")

t2 = time()
println("Time to load: ", t2 - t1)

fcode = pyimport("fcode")
t3 = time()
println("Time to load fcode: ", t3 - t2)
cfg = fcode.config

prefix = cfg.prefix
np_BZ = np.array(cfg.brillouin_zone)
projections = cfg.projections

nx, ny, nz = cfg.k_mesh
dim = cfg.dimension
println("Dimension: ", dim)
if dim == 2
    nz = 1
end
nk = nx * ny * nz
nw = cfg.w_pts 
U = cfg.onsite_U
BZ = cfg.brillouin_zone
mu = cfg.fermi_energy
t4 = time()
println("Time to load config: ", t4 - t3)

const beta = 1/cfg.Temperature
println("Beta: ", beta)
const pi = π


wD = 0.2

function to_IBZ(k)
    tolerance = 1e-10  # small tolerance to account for floating-point errors
    q = abs.(k)
    q .= ifelse.(abs.(q .- 2π) .< tolerance, 0.0, ifelse.(q .> π, -(q .- 2π), q))
    q = abs.(q)
    q .+= 1e-4
    q .= ifelse.(q .> π, q .- 1e-4, q)
    return q
end

function epsilon(k)
    if dim == 2
        #return k[1]^2 + k[2]^2
        return -2 * (cos(k[1]) + cos(k[2]))
    end
    return -2 * (cos(k[1]) + cos(k[2]) + cos(k[3]))
end

function V(k1, k2, w)
    qm = k1 - k2
    qm = to_IBZ(qm)
    qm .= round.(qm, digits=6)
    iw = round(imag(w), digits=6)
    Xm = Susceptibility(qm, iw)
    Vm = U^2 * Xm / (1 - U*Xm) + U^3 * Xm^2 / (1 - U^2 * Xm^2)
    return Vm
end

function get_bandwidth()
    maxval = -1000
    minval = 1000
    for i in 1:200, j in 1:200, k in 1:200
        kvec = get_kvec(i, j, k)
        eps = epsilon(kvec)
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
    e = epsilon(k) - mu
    lambda = -1.0
    if abs(e) < wD
        return lambda 
    end
    return 0.0
end

function fill_V_arr(iw, frequency_dependence)
    Vw_arr = Array{ComplexF64}(undef, length(iw), nx, ny, nz)
    @threads for i in 1:nx
        @inbounds for j in 1:ny 
            @inbounds for k in 1:nz
                @inbounds for l in 1:length(iw)
                    kvec = get_kvec(i, j, k)
                    w1 = iw[l]
                    #Vw_arr[l, i, j] = V(kvec, zeros(dim), w1)
                    #Vw_arr[l, i, j, k] = paper2_V(w1)
                    Vw_arr[l, i, j, k] = BCS_V(kvec, w1)
                end
            end
        end
    end
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

function get_kvec(i, j, k)
    temp = BZ * [i / nx - 0.5, j / ny - 0.5, k / nz - 0.5]
    return temp[1:dim]
end

function condense_to_F_and_G!(phi, Z, chi, F_arr, G_arr, iw, sigma, mu, projection)
    @threads for l in 1:length(iw)
        @inbounds for i in 1:nx
            @inbounds for j in 1:ny
                @inbounds for k in 1:nz
                    w = iw[l]
                    kvec = get_kvec(i, j, k)
                    phi_el = phi[l] * projection[i, j, k]
                    denom = get_denominator(imag(w), phi_el, Z[l, i, j, k], chi[l, i, j, k] + epsilon(kvec) - mu)
                    F_arr[l, i, j, k] = -phi_el / denom
                    G_arr[l, i, j, k] = -(w * Z[l, i, j, k] + epsilon(kvec) - mu + chi[l, i, j, k]) / denom
                end
            end
        end
    end
end

function update!(F, G, V_arr, phi, Z, chi, iw, sigma, projection)
    F_rt = kw_to_rtau(F, 'F', mesh)
    G_rt = kw_to_rtau(G, 'F', mesh)

    phit = V_arr .* F_rt / 2
    sigmat = -V_arr .* G_rt / 2

    phi_full = rtau_to_kw(phit, 'F', mesh)
    sigma .= rtau_to_kw(sigmat, 'F', mesh)

    phi .= phi_full[:, 1, 1, 1] / projection[1, 1, 1]
    fill_Z_chi(iw, sigma, Z, chi)
end

function update_fourier_transform!(F, G, V_arr, phi, Z, chi, iw, sigma, projection)
    F_rt = fft(F)
    G_rt = fft(G)

    phit = V_arr .* F_rt 
    sigmat = -V_arr .* G_rt 

    phi_full = fftshift(ifft(phit)) / (nw * nx * ny * nz)
    sigma .= fftshift(ifft(sigmat)) / (nw * nx * ny * nz)

    phi .= phi_full[:, 1, 1, 1] ./ projection[1, 1, 1]
    fill_Z_chi(iw, sigma, Z, chi)
end

function fill_Z_chi(iw, sigma, Z, chi)
    Z .= 1.0
    chi .= 0.0
    return
    @threads for i in 1:nw
        @inbounds for j in 1:nx
            @inbounds for k in 1:ny
                @inbounds for l in 1:nz
                    sp, sm = sigma[i, j, k, l], sigma[nw - i + 1, j, k, l]
                    Z[i, j, k, l] = 1.0 - 0.5 * (sp - sm) / iw[i]
                    chi[i, j, k, l] = 0.5 * (sp + sm)
                end
            end
        end
    end
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
    IR_tol = 1e-10
    scf_tol = 1e-6

    dosfile = prefix * "_DOS.dat"
    DOS = Field.Field_R(dosfile)
    println("fermi energy: ", mu)
    n0 = DOS(mu)
    println("DOS: ", n0)
    #global Susceptibility = CondensedMatterField.CMF(input_data_file)
    println("Loaded Susceptibility")

    frequency_dependence = true
    #if V([0.2,0.3,0.4], 0) == V([0.2,0.3,0.4], 0.3)
    #    frequency_dependence = false
    #end
    frequency_dependence = false
    if frequency_dependence
        println("Creating IRMesh")
        global mesh = IR_Mesh(IR_tol)
        global fnw = mesh.fnw
        global fntau = mesh.fntau
        global bnw = mesh.bnw
        global bntau = mesh.bntau
        println("IRMesh created")
        iw = Array{ComplexF64}(undef, fnw)
        iv = Array{ComplexF64}(undef, bnw)
        for i in 1:fnw iw[i] = valueim(mesh.IR_basis_set.smpl_wn_f.sampling_points[i], beta) end
        for i in 1:bnw iv[i] = valueim(mesh.IR_basis_set.smpl_wn_b.sampling_points[i], beta) end

        println("Filled iw and iv")
        println("Min & Max iw: ", minimum(imag.(iw)), " ", maximum(imag.(iw)))
        global nw = length(iw)
    else
        iw = Array{ComplexF64}(undef, nw)
        iv = Array{ComplexF64}(undef, nw)
        for i in 1:nw
            n = -nw / 2 + i - 1
            iw[i] = Complex(0.0, (2 * n - 1) * pi / beta)
            iv[i] = Complex(0.0, 2 * n * pi / beta)
        end
    end

    println("Initializing projection")
    projection = zeros(Float64, nx, ny, nz)
    if projections == "s"
        projection .= 1.0
    elseif projections == "d"
        for i in 1:nx, j in 1:ny, k in 1:nz
            kvec = get_kvec(i, j, k)
            projection[i, j, k] = cos(kvec[1]) - cos(kvec[2])
        end
    end

    println("Initializing phi, Z, and chi")
    phi_arr = ComplexF64.(0.001 ./ (imag.(iw).^2 .+ 1))
    Z_arr = ones(Complex{Float64}, nw, nx, ny, nz)
    chi_arr = zeros(Complex{Float64}, nw, nx, ny, nz)


    println("Initializing V")
    sigma = zeros(Complex{Float64}, nw, nx, ny, nz)
    V_arr = fill_V_arr(iv, frequency_dependence)

    println("Initializing F and G")
    F_arr = Array{ComplexF64}(undef, nw, nx, ny, nz)
    G_arr = Array{ComplexF64}(undef, nw, nx, ny, nz)
    condense_to_F_and_G!(phi_arr, Z_arr, chi_arr, F_arr, G_arr, iw, sigma, mu, projection)

    Ne = round(get_number_of_electrons(G_arr), digits=6)

    phierr = 0.0
    prev_phi_arr = copy(phi_arr)
    println("Starting Convergence Loop")
    iterations = 1000
    #iterations = 1
    for i in 1:iterations
        if frequency_dependence
            update!(F_arr, G_arr, V_arr, phi_arr, Z_arr, chi_arr, iw, sigma, projection)
        else
            update_fourier_transform!(F_arr, G_arr, V_arr, phi_arr, Z_arr, chi_arr, iw, sigma, projection)
        end

        #phierr = round(sum(abs.(prev_phi_arr - phi_arr)) / (nw * nx * ny * nz),digits=8)
        phierr = maximum(abs.(phi_arr - prev_phi_arr))
        max_phi = round(maximum(abs.(phi_arr)), digits=6)
        #print("Iteration $i: Max Phi = $max_phi, mu = $mu, Error = $phierr           \r")
        print("Iteration $i: Error = $phierr           \r")

        if abs(phierr) < scf_tol || max_phi < 1e-10
            println("Converged after ", i, " iterations                                                                   ")
            break
        end
        prev_phi_arr = copy(phi_arr)
        global mu = find_chemical_potential(mu, phi_arr, Z_arr, chi_arr, F_arr, G_arr, iw, sigma, Ne, projection)
        condense_to_F_and_G!(phi_arr, Z_arr, chi_arr, F_arr, G_arr, iw, sigma, mu, projection)
    end

    max_phi = maximum(abs.(phi_arr))
    max_Z = maximum(abs.(Z_arr))

    println("Error: ", phierr, "            ")
    println("Max phi: ", max_phi)
    println("Max Z: ", max_Z)
    lambda = 0.2
    println("Analytical Phi max: ", wD / sinh(( n0 / lambda)))
    reordered_phi = Array{ComplexF64}(undef, nx, ny, nz, nw)
    reordered_Z = Array{ComplexF64}(undef, nx, ny, nz, nw)
    points = Matrix{Float64}(undef, nx*ny*nz*nw, dim+1)
    @threads for i in 1:nw
        for j in 1:nx, k in 1:ny, l in 1:nz
            reordered_phi[j, k, l, i] = phi_arr[i] * projection[j, k, l]
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

    p = plot(npts, real.(phi_arr), label="Phi", xlims=(-500, 500))
    #plot!(p, phi_average_imag, label="Imag Phi")
    display(p)
    #readline()
    p = plot(npts, real.(Z_w), label="Z", xlims=(-500, 500))
    display(p)
    #readline()

    println("Saving phi and Z")
    #save_data(prefix*"_phi.dat", points, reordered_phi, dim, true, true, false)
    #save_data(prefix*"_Z.dat", points, reordered_Z, dim, true, true, false)
end 

function get_denominator(w1, phi_el, Z_el, eps)
    return abs2(w1 * real(Z_el)) + abs2(eps) + abs2(real(phi_el))
end

function average_over_k(arr)
    #return Float64.(arr[1, 1, 1, :])
    average = zeros(Float64, length(arr[1, 1, 1, :]))
    for i in nx, j in ny, k in nz
        average .+= (arr[i, j, k, :])
    end
    return average
end

@inline function F_gpu(iw, e, phi) 
    return phi / (-iw^2 + phi^2 + e^2)
end
@inline function V_gpu(iv, e, wD, lambda) 
    return abs(e) > wD ? 0.0 : lambda
    #return lambda
    return lambda * 2 * wD^2 / (wD^2 - iv^2)
end
@inline function epsilon_cuda(kx, ky, kz) 
    return kx^2 + ky^2
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

    lambda = 1.0
    wD = 0.2
    mu = 1.0

    F(iw, e, phi) = phi / (-iw^2 + phi^2 + e^2)
    V(iv, e) = abs(e) > wD ? 0.0 : lambda

    F_arr = CUDA.zeros(ComplexF64, pts, nx, ny)
    V_arr = CUDA.zeros(ComplexF64, pts, nx, ny)
    e_arr = CUDA.zeros(Float64, nx, ny)
    w_arr = CUDA.zeros(ComplexF64, pts)
    phi = CUDA.ones(Float64, pts) * 1e-3
    BZ_GPU = CuArray(BZ)
    println("Filling arrays")
    threads = 256
    blocks_per_grid = cld(pts, threads)
    CUDA.@cuda threads=threads blocks=blocks_per_grid fill_gpu!(F_arr, V_arr, e_arr, phi, pts, nx, ny, beta, BZ_GPU, mu, w_arr, F, V, epsilon_cuda)
    
    println("Starting Fourier Transforms")
    error = CUDA.ones(Float64, pts)
    max_arr = CUDA.zeros(Float64, pts)
    for i in 1:50
        F_rt = CUDA.CUFFT.fft(F_arr)
        V_rt = CUDA.CUFFT.fft(V_arr)
        phi_rt = CUDA.dot(F_rt, V_rt)
        ifft = CUDA.CUFFT.ifft(phi_rt)
        result = CUDA.fftshift(CUDA.CUFFT.ifft(phi_rt)) / (pts * nk)

        CUDA.@cuda threads=threads blocks=blocks_per_grid error_gpu(result, phi, error, max_arr, pts, nx, ny)
        error_total = sum(Array(error))
        max_total = max(Array(max_arr))
        println("Error total: ", error_total)
        println("Max val: ", max_total)

        if error_total < 1e-4
            break
        end
        CUDA.@cuda threads=threads blocks=blocks_per_grid update_gpu!(F_arr, e_arr, result, pts, nx, ny, w_arr, F)
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
    result = zeros(ComplexF64, pts)
    F_arr = zeros(ComplexF64, pts)
    V_arr = zeros(ComplexF64, pts)

    w_arr = Array{Float64}(undef, pts)
    iw_arr = Array{Float64}(undef, fnw)
    for i in 1:fnw
        iw_arr[i] = imag(valueim(basis.IR_basis_set.smpl_wn_f.sampling_points[i], beta))
    end

    F_ir = Array{ComplexF64}(undef, fnw)
    V_ir = Array{ComplexF64}(undef, bnw)
    phi = 1.0
    lambda = 1.0
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
            result[i] += ( V(iwn - iwm) * F(iwm) ) / beta
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
    result = zeros(ComplexF64, pts)
    F_arr = zeros(ComplexF64, pts, nx, ny, nz)
    V_arr = zeros(ComplexF64, pts, nx, ny, nz)

    w_arr = Array{Float64}(undef, pts)
    iw_arr = Array{Float64}(undef, fnw)
    for i in 1:fnw
        iw_arr[i] = imag(valueim(basis.IR_basis_set.smpl_wn_f.sampling_points[i], beta))
    end

    F_ir = Array{ComplexF64}(undef, fnw, nx, ny, nz)
    V_ir = Array{ComplexF64}(undef, bnw, nx, ny, nz)
    phi = 1.0
    lambda = 1.0
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
                        result[i] += ( V(iwn - iwm) * F(iwm, em) ) / (beta * nk)
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


@inline function F_gpu(iw, e, phi)
    return phi / (-iw^2 + phi^2 + e^2)
end
@inline function V_gpu(iv, e, wD, lambda)
    return abs(e) > wD ? 0.0 : lambda
    #return lambda
    return lambda * 2 * wD^2 / (wD^2 - iv^2)
end
@inline function epsilon_cuda(kx, ky, kz)
    return -2 * (cos(kx) + cos(ky))
    return kx^2 + ky^2
end

function compute_gpu_sum!(result, pts, nx, ny, nz, beta, phi, wD, lambda, w_arr_CUDA, BZ, mu)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x

    if i ≤ pts
        iwn = ComplexF64(w_arr_CUDA[i])
        sum = 0

        for j in 1:pts, a in 1:nx, b in 1:ny, c in 1:nz
            iwm = ComplexF64(w_arr_CUDA[j])
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

function VF_4sum_GPU(phi, pts, phi_ir, basis)
    println("Starting")
    fnw, bnw, fntau, bntau = basis.fnw, basis.bnw, basis.fntau, basis.bntau
    nk = nx * ny * nz

    lambda = 1.0
    wD = 0.2
    mu = 1.0

    F(iw, e, phi) = phi / (-iw^2 + phi^2 + e^2)
    #V(iv) = lambda * 2 * wD^2 / (wD^2 - iv^2)
    V(iv, e) = abs(e) > wD ? 0.0 : lambda

    result = CUDA.zeros(ComplexF64, pts)
    F_arr = zeros(ComplexF64, pts, nx, ny, nz)
    V_arr = zeros(ComplexF64, pts, nx, ny, nz)
    w_arr = zeros(ComplexF64, pts)
    w_array = Array{ComplexF64}(undef, pts)
    println("Filling arrays")
    @threads for i in 1:pts
        for j in 1:nx, k in 1:ny, l in 1:nz
            n = -pts / 2 + i - 1
            iwn = Complex(0.0, (2 * n + 1) * pi / beta)
            iwm = Complex(0.0, (2 * n) * pi / beta)
            w_arr[i] = iwn
            w_array[i] = iwn
            kx = BZ[1, 1] * (j / nx - 0.5) + BZ[1, 2] * (k / ny - 0.5) + BZ[1, 3] * (l / nz - 0.5)
            ky = BZ[2, 1] * (j / nx - 0.5) + BZ[2, 2] * (k / ny - 0.5) + BZ[2, 3] * (l / nz - 0.5)
            kz = BZ[3, 1] * (j / nx - 0.5) + BZ[3, 2] * (k / ny - 0.5) + BZ[3, 3] * (l / nz - 0.5)
            e = epsilon_cuda(kx, ky, kz) - mu
            e = epsilon(get_kvec(j, k, l)) - mu
            F_arr[i, j, k, l] = F(iwn, e, phi[i])
            V_arr[i, j, k, l] = V(iwm, e)
        end
    end
    w_arr_CUDA = CuArray(w_array)
    BZ_CUDA = CuArray(BZ)

    iw_arr = Array{Float64}(undef, fnw)
    for i in 1:fnw
        iw_arr[i] = imag(valueim(basis.IR_basis_set.smpl_wn_f.sampling_points[i], beta))
    end

    F_ir = Array{ComplexF64}(undef, fnw, nx, ny, nz)
    V_ir = Array{ComplexF64}(undef, bnw, nx, ny, nz)

    println("Filling IR arrays")
    for i in 1:fnw, j in 1:nx, k in 1:ny, l in 1:nz
        iw = valueim(basis.IR_basis_set.smpl_wn_f.sampling_points[i], beta)
        F_ir[i, j, k, l] = F(iw, epsilon(get_kvec(j, k, l)) - mu, phi_ir[i])
    end
    for i in 1:bnw, j in 1:nx, k in 1:ny, l in 1:nz
        iv = valueim(basis.IR_basis_set.smpl_wn_b.sampling_points[i], beta)
        V_ir[i, j, k, l] = V(iv, epsilon(get_kvec(j, k, l)) - mu)
    end

    println("Starting GPU sum")
    phi_CUDA = CuArray(phi)

    CUDA.@sync @cuda threads=512 blocks=cld(pts, 512) compute_gpu_sum!(result, pts, nx, ny, nz, beta, phi_CUDA, wD, lambda, w_arr_CUDA, BZ_CUDA, mu)
#    CUDA.@sync @cuda threads=threads_per_block blocks=blocks_per_grid compute_gpu_sum!(
#                                                                                       result, pts, nx, ny, nz, beta, phi, wD, lambda, w_arr_CUDA, BZ_CUDA)


    cpu_result = Array(result)

    println("Starting Fourier Transforms")
    F_rt = fft(F_arr)
    V_rt = fft(V_arr)
    phi_rt = F_rt .* V_rt
    result2 = fftshift(ifft(phi_rt)) / (pts * nk)

    F_irt = kw_to_rtau(F_ir, 'F', basis)
    V_irt = kw_to_rtau(V_ir, 'B', basis)
    phi_irt = F_irt .* V_irt
    result3 = rtau_to_kw(phi_irt, 'F', basis)

    println("Printing maxvals")
    println("Max result: ", maximum(abs.(cpu_result)))
    println("Max result2: ", maximum(abs.(result2)))
    println("Max result3: ", maximum(abs.(result3)))

    println("Plotting")
    result_plot = cpu_result
    max2_index = argmax(abs.(result2))
    println("Max2 index: ", max2_index)
    max2_indices = Tuple(CartesianIndices(result2)[max2_index])
    result2_plot = result2[:, max2_indices[2], max2_indices[3], max2_indices[4]]
    result3_plot = sum(result3, dims=(2, 3, 4))[:, 1, 1, 1] / (nx * ny * nz)

    #p = plot(imag.(w_arr), real.(result_plot), label="Summed", xlims=(-20, 20))
    #plot!(p, imag.(w_arr), real.(result2_plot), label="Ft'd", xlims=(-20, 20))
    ##plot!(p, iw_arr, real.(result3_plot), label="IR'd", xlims=(-20, 20))
    #display(p)
    #readline()
    return cpu_result, result3_plot
end

function gpu_eliashberg()
    pts = 10000
    phi = ones(ComplexF64, pts) * 1e-3
    basis = IR_Mesh(1e-10)
    fnw, bnw, fntau, bntau = basis.fnw, basis.bnw, basis.fntau, basis.bntau
    phi_ir = ones(ComplexF64, fnw)
    prev_phi = ones(ComplexF64, pts)
    for i in 1:20
        prev_phi .= phi
        println("Iteration: ", i)
        phi, phi_ir = Eliashberg.VF_4sum_GPU(phi, pts, phi_ir, basis)
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

end # module

#Eliashberg.test_max_bcs_no_w()
Eliashberg.eliashberg_global()
#Eliashberg.gpu_eliashberg()
#Eliashberg.VF_3sum_GPU(Int(1e5))
