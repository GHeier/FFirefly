module Eliashberg

t1 = time()
include("../objects/mesh.jl")
using .IRMesh

using CUDA, FFTW
using Plots
using Roots
using SparseIR
import SparseIR: Statistics, value, valueim
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


wD = 0.1

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
    return -1*lambda * v^2 / (n^2 + v^2)
end

function BCS_V(k, w)
    e = epsilon(k) - mu
    lambda = -0.2
    if abs(e) < wD
        return lambda / (1 + 0.01 * imag(w)^2)
    end
    return lambda / (1 + 0.01 * imag(w)^2) / 100
end

function fill_V_arr!(V_arr, iw)
    Vw_arr = Array{ComplexF64}(undef, bnw, nx, ny, nz)
    @threads for i in 1:nx
        @inbounds for j in 1:ny 
            @inbounds for k in 1:nz
                @inbounds for l in 1:bnw
                    kvec = get_kvec(i, j, k)
                    w1 = iw[l]
                    #Vw_arr[l, i, j] = V(kvec, zeros(dim), w1)
                    #Vw_arr[l, i, j, k] = paper2_V(w1)
                    Vw_arr[l, i, j, k] = BCS_V(kvec, w1)
                end
            end
        end
    end
    V_arr .= kw_to_rtau(Vw_arr, 'B', mesh)
    println("Filled V_arr")
end

function initialize_phi_Z_chi!(phi_arr, Z_arr, chi_arr, iw)
    @threads for i in 1:fnw
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

function condense_to_F_and_G!(phi, Z, chi, F_arr, G_arr, iw, sigma, mu)
    @threads for l in 1:fnw
        @inbounds for i in 1:nx
            @inbounds for j in 1:ny
                @inbounds for k in 1:nz
                    w = iw[l]
                    kvec = get_kvec(i, j, k)
                    denom = get_denominator(imag(w), phi[l, i, j, k], Z[l, i, j, k], chi[l, i, j, k] + epsilon(kvec) - mu)
                    F_arr[l, i, j, k] = -phi[l, i, j, k] / denom
                    G_arr[l, i, j, k] = -(w * Z[l, i, j, k] + epsilon(kvec) - mu + chi[l, i, j, k]) / denom
                end
            end
        end
    end
end

function update!(F, G, V_arr, phi, Z, chi, iw, sigma)
    F_rt = kw_to_rtau(F, 'F', mesh)
    G_rt = kw_to_rtau(G, 'F', mesh)

    phit = V_arr .* F_rt 
    sigmat = -V_arr .* G_rt 

    phi .= rtau_to_kw(phit, 'F', mesh)
    sigma .= rtau_to_kw(sigmat, 'F', mesh)

    @threads for i in 1:fnw
        @inbounds for j in 1:nx
            @inbounds for k in 1:ny
                @inbounds for l in 1:nz
                    w = iw[i]
                    sp, sm = sigma[i, j, k, l], sigma[fnw - i + 1, j, k, l]
                    Z[i, j, k, l] = 1.0 - 0.5 * (sp - sm) / w
                    #Z[i, j, k, l] = 1
                    chi[i, j, k, l] = 0.5 * (sp + sm)
                end
            end
        end
    end
end


function find_chemical_potential(mu_initial, phi_arr, Z_arr, chi_arr, F_arr, G_arr, iw, sigma, target_Ne)
    return mu_initial
    f(mu) = begin
        condense_to_F_and_G!(phi_arr, Z_arr, chi_arr, F_arr, G_arr, iw, sigma, mu)
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

function eliashberg_sparse_ir()
    println("Beginning Eliashberg")
    IR_tol = 1e-10
    scf_tol = 1e-4
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
    input_data_file = prefix * "_ckio_ir.dat"
    println("input_data_file: ", input_data_file)
    #global Susceptibility = CondensedMatterField.CMF(input_data_file)
    println("Loaded Susceptibility")

    phi_arr = Array{ComplexF64}(undef, fnw, nx, ny, nz)
    Z_arr = Array{ComplexF64}(undef, fnw, nx, ny, nz)
    chi_arr = Array{ComplexF64}(undef, fnw, nx, ny, nz)
    println("Initializing phi, Z, and chi")
    initialize_phi_Z_chi!(phi_arr, Z_arr, chi_arr, iw)

    println("Initializing V")
    sigma = zeros(Complex{Float64}, fnw, nx, ny, nz)
    V_arr = Array{ComplexF64}(undef, bntau, nx, ny, nz)
    fill_V_arr!(V_arr, iv)

    println("Initializing F and G")
    F_arr = Array{ComplexF64}(undef, fnw, nx, ny, nz)
    G_arr = Array{ComplexF64}(undef, fnw, nx, ny, nz)
    condense_to_F_and_G!(phi_arr, Z_arr, chi_arr, F_arr, G_arr, iw, sigma, mu)

    Ne = round(get_number_of_electrons(G_arr), digits=6)

    phierr = 0.0
    prev_phi_arr = copy(phi_arr)
    println("Starting Convergence Loop")
    #iterations = 1000
    iterations = 1
    for i in 1:iterations
        update!(F_arr, G_arr, V_arr, phi_arr, Z_arr, chi_arr, iw, sigma)

        phierr = round(sum(abs.(prev_phi_arr - phi_arr)) / (fnw * nx * ny * nz),digits=8)
        max_phi = round(maximum(abs.(phi_arr)), digits=6)
        #print("Iteration $i: Max Phi = $max_phi, mu = $mu, Error = $phierr           \r")
        print("Iteration $i: Error = $phierr           \r")

        if abs(phierr) < scf_tol || max_phi < 1e-10
            println("Converged after ", i, " iterations                                                                   ")
            break
        end
        prev_phi_arr = copy(phi_arr)
        global mu = find_chemical_potential(mu, phi_arr, Z_arr, chi_arr, F_arr, G_arr, iw, sigma, Ne)
        condense_to_F_and_G!(phi_arr, Z_arr, chi_arr, F_arr, G_arr, iw, sigma, mu)
    end

    max_phi = maximum(abs.(phi_arr))
    max_Z = maximum(abs.(Z_arr))

    println("Error: ", phierr, "            ")
    println("Max phi: ", max_phi)
    println("Max Z: ", max_Z)
    reordered_phi = Array{ComplexF64}(undef, nx, ny, nz, fnw)
    reordered_Z = Array{ComplexF64}(undef, nx, ny, nz, fnw)
    points = Matrix{Float64}(undef, nx*ny*nz*fnw, dim+1)
    for i in 1:fnw, j in 1:nx, k in 1:ny, l in 1:nz
        reordered_phi[j, k, l, i] = phi_arr[i, j, k, l]
        reordered_Z[j, k, l, i] = Z_arr[i, j, k, l]
        w = iw[i]
        kvec = get_kvec(j, k, l)
        idx = (i - 1) * nx * ny * nz + (j - 1) * ny * nz + (k - 1) * nz + l
        points[idx, 1:dim] .= kvec
        points[idx, dim+1] = imag(w)
    end

    npts = Array{Integer}(undef, fnw)
    for i in 1:fnw
        w = iw[i]
        n = trunc((beta * imag(w) / pi - 1) / 2)
        npts[i] = n
    end
    phi_average = average_over_k(real.(reordered_phi))
    Z_average = average_over_k(real.(reordered_Z))
    phi_average_imag = average_over_k(imag.(reordered_phi))

    p = plot(npts, phi_average, label="Phi", xlims=(-500, 500))
    #plot!(p, phi_average_imag, label="Imag Phi")
    display(p)
    readline()
    p = plot(npts, Z_average, label="Z", xlims=(-500, 500))
    display(p)
    readline()

    println("Saving phi and Z")
    #save_arr_as_meshgrid(points, reordered_phi, "phi.dat", true)
    #save_arr_as_meshgrid(points, reordered_Z, "Z.dat", true)
end 

function average_over_k(arr)
    #return Float64.(arr[1, 1, 1, :])
    average = zeros(Float64, fnw)
    for i in nx, j in ny, k in nz
        average .+= (arr[i, j, k, :])
    end
    return average
end

function eliashberg_isotropic()
    println("Beginning Eliashberg(projection)")
    IR_tol = 1e-10
    scf_tol = 1e-4
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
    input_data_file = prefix * "_ckio_ir.dat"
    println("input_data_file: ", input_data_file)
    #global Susceptibility = CondensedMatterField.CMF(input_data_file)
    println("Loaded Susceptibility")

    phi_arr = Array{ComplexF64}(undef, fnw)
    phi_arr .= 0.001 / (imag.(iw).^2 + 1)
    Z_arr = ones(Complex{Float64}, fnw)
    chi_arr = zeros(Complex{Float64}, fnw)

    sigma = zeros(Complex{Float64}, fnw, nx, ny, nz)
    V_arr = Array{ComplexF64}(undef, bntau, nx, ny, nz)
    fill_V_arr!(V_arr, iv)

    F_arr = Array{ComplexF64}(undef, fnw, nx, ny, nz)
    G_arr = Array{ComplexF64}(undef, fnw, nx, ny, nz)
    condense_to_F_and_G!(phi_arr, Z_arr, chi_arr, F_arr, G_arr, iw, sigma, mu)
    
    Ne = round(get_number_of_electrons(G_arr), digits=6)

    phierr = 0.0
    prev_phi_arr = copy(phi_arr)
    println("Starting Convergence Loop")
    iterations = 1000
    for i in 1:iterations
        update_isotropic!(F_arr, G_arr, V_arr, phi_arr, Z_arr, chi_arr, iw, sigma)

        phierr = round(sum(abs.(prev_phi_arr - phi_arr)) / fnw, digits=8)
        max_phi = round(maximum(abs.(phi_arr)), digits=6)
        print("Iteration $i: Max Phi = $max_phi, mu = $mu, Error = $phierr           \r")

        if abs(phierr) < scf_tol || max_phi < 1e-10
            println("Converged after ", i, " iterations                                                                   ")
            break
        end
        prev_phi_arr = copy(phi_arr)
        global mu = find_chemical_potential(mu, phi_arr, Z_arr, chi_arr, F_arr, G_arr, iw, sigma, Ne)
        condense_to_F_and_G!(phi_arr, Z_arr, chi_arr, F_arr, G_arr, iw, sigma, mu)
    end

end

function update_isotropic!(F, G, V_arr, phi, Z, chi, iw, sigma)
    F_rt = kw_to_rtau(F, 'F', mesh)
    G_rt = kw_to_rtau(G, 'F', mesh)

    phit = Array{ComplexF64}(undef, fntau, nx, ny, nz)
    sigmat = Array{ComplexF64}(undef, fntau, nx, ny, nz)
    for i in 1:fntau, j in 1:nx, k in 1:ny, l in 1:nz
        phit[i, j, k, l] = V_arr[i, j, k, l] * F_rt[i, j, k, l]
        sigmat[i, j, k, l] = -V_arr[i, j, k, l] * G_rt[i, j, k, l]
    end

    phi .= rtau_to_kw(phit, 'F', mesh)
    sigma .= rtau_to_kw(sigmat, 'F', mesh)

    @threads for i in 1:fnw
        @inbounds for j in 1:nx
            @inbounds for k in 1:ny
                @inbounds for l in 1:nz
                    w = iw[i]
                    sp, sm = sigma[i, j, k, l], sigma[fnw - i + 1, j, k, l]
                    Z[i] = 1.0 - 0.5 * (sp - sm) / w
                    chi[i] = 0.5 * (sp + sm)
                end
            end
        end
    end
end

function condense_to_F_and_G!(phi_arr, Z_arr, chi_arr, F_arr, G_arr, iw, sigma_arr, mu, projections)
    @threads for l in 1:fnw
        @inbounds for i in 1:nx
            @inbounds for j in 1:ny
                @inbounds for k in 1:nz
                    w = iw[l]
                    kvec = get_kvec(i, j, k)
                    phi = phi_arr[l] * projections[l]
                    Z = Z_arr[l]
                    chi = chi_arr[l]
                    denom = get_denominator(imag(w), phi, Z, chi + epsilon(kvec) - mu)
                    F_arr[l, i, j, k] = -phi[l] / denom
                    G_arr[l, i, j, k] = -(w * Z[l] + epsilon(kvec) - mu + chi[l]) / denom
                end
            end
        end
    end
end


function get_denominator(w1, phi_el, Z_el, eps)
    #return (abs2(w1 * Z_el) + abs2(phi_el))^(0.5) / pi
    #return (abs(w1 * Z_el) + abs(phi_el)) 
    return abs2(w1 * real(Z_el)) + abs2(eps) + abs2(real(phi_el))
end

function gpu_sum()
    # Initialize data on GPU
    scf_tol = 1e-4
    fnw = 1000
    npts = Array{Integer}(undef, fnw)
    iw = Array{ComplexF64}(undef, fnw)
    phi = Array{ComplexF64}(undef, fnw)
    Z = Array{ComplexF64}(undef, fnw)
    for i in 1:fnw
        #npts[i] = -fnw/2 + i - 1
        npts[i] = i
        iw[i] = Complex(0.0, (pi*(2*(npts[i] - 1) + 1) / beta))
        phi[i] = 1 / abs(iw[i].im^2 + 1)
        Z[i] = 1 / abs(iw[i].im^2 + 1)
    end

    iw_gpu = CuArray(iw)
    phi_gpu = CuArray(phi)
    Z_gpu = CuArray(Z)

    BZ_gpu = CuArray(BZ)

    for iters in 1:20
        print("Iteration $iters: ")
        new_phi = CuArray(zeros(Complex{Float64}, fnw))
        new_Z = CuArray(zeros(Complex{Float64}, fnw))

        #@device_code_warntype kernel_summation!(
        #    new_phi, new_Z, phi_gpu, Z_gpu, iw_gpu, nx, ny, beta, BZ_gpu, mu
        #)
        # GPU Kernel for summation
        @cuda threads=512 blocks=cld(fnw, 512) kernel_summation!(
             new_phi, new_Z, phi_gpu, Z_gpu, iw_gpu, nx, ny, nz, beta, BZ_gpu, mu
        )

        # Convergence check
        error = CUDA.maximum(abs.(new_phi - phi_gpu))
        println("Error = ", error)
        if error < scf_tol
            println("Converged after $iters iterations")
            break
        end

        # Update variables
        phi_gpu .= new_phi
        Z_gpu .= new_Z
    end
    phi = Array(phi_gpu)
    Z = Array(Z_gpu)
    println("Max phi: ", maximum(abs.(phi)))
    println("Max Z: ", maximum(abs.(Z)))
    p = plot(npts, real.(phi))
    display(p)
    readline()
    p = plot(npts, real.(Z))
    display(p)
    readline()
end

# GPU Kernel Definition
function kernel_summation!(
    new_phi::CuDeviceVector{ComplexF64},
    new_Z::CuDeviceVector{ComplexF64},
    phi::CuDeviceVector{ComplexF64},
    Z::CuDeviceVector{ComplexF64},
    iw::CuDeviceVector{ComplexF64},
    nx::Int,
    ny::Int,
    nz::Int,
    beta::Float64,
    BZ::CuDeviceMatrix{Float64},
    mu::Float64
)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    if i > length(new_phi)
        return
    end

    temp_phi = ComplexF64(0.0, 0.0)
    temp_Z = ComplexF64(1.0, 0.0)

    #V(x) = 1.0 / (imag(x)^2 + 0.01)
    V(x) = paper2_V(x)
    e(kx, ky, kz) = -2 * (cos(kx) + cos(ky) + cos(kz))
    get_denom(iw, phi, Z, eps) = abs2(imag(iw)*Z) + abs2(eps) + abs2(phi)

    for j in 1:length(phi)
        for k in 1:nx, l in 1:ny, m in 1:nz
            kx = BZ[1, 1] * (k / nx) + BZ[1, 2] * (l / ny) + BZ[1, 3] * (m / nz) 
            ky = BZ[2, 1] * (k / nx) + BZ[2, 2] * (l / ny) + BZ[2, 3] * (m / nz) 
            kz = BZ[3, 1] * (k / nx) + BZ[3, 2] * (l / ny) + BZ[3, 3] * (m / nz) 
            denom = get_denom(iw[j], phi[j], Z[j], e(kx, ky, kz) - mu)
            w = iw[i] - iw[j]
            w2 = iw[i] + iw[j] + 1.0 * pi / beta

            V_val = V(w)
            V_val2 = V(w2)
            # Accumulate results
            temp_phi -= (1 / beta) * (V_val + V_val2) * phi[j] / denom
            temp_Z -= (1 / beta) * (V_val - V_val2) * (iw[j] / iw[i]) * Z[j] / denom
        end
    end

    new_phi[i] = temp_phi / (nx * ny * nz)
    new_Z[i] = temp_Z / (nx * ny * nz)
    return nothing
end

function phonon_V(w)
    # Define constants
    g = 1.0
    wq = 1.0
    e = 0.4

    # Define functions
    V = -g^2 * 2 * wq / (-wq^2 + w^2)
    return V
end

function Greens_function(w)
    e = 0.4
    return 1 / (w - e)
end

function test_4d_multiply()
    IR_basis_set = FiniteTempBasisSet(beta, 100.0, 1e-15)
    mesh = IRMesh.Mesh(IR_basis_set, nx, ny, nz)
    fnw, bnw = mesh.fnw, mesh.bnw

    # Define constants
    g = 1.0
    wq = 1.0
    e = 0.4

    # Define functions
    V(w) = -g^2 * 2 * wq / (-wq^2 + w^2)
    G(w) = 1 / (w - e)
    n(e) = 1 / (exp(beta * e) - 1)
    f(e) = 1 / (exp(beta * e) + 1)
    Sigma(w) = g^2 * ( (1 + n(wq) - f(e)) / (w - (e + wq)) + (n(wq) + f(e)) / (w - (e - wq)))

    # Allocate arrays for storing functions
    Vw = Array{ComplexF64}(undef, bnw, nx, ny, nz)
    Gw = Array{ComplexF64}(undef, fnw, nx, ny, nz)
    Sigmaw = Array{ComplexF64}(undef, fnw, nx, ny, nz)
    iw = Array{ComplexF64}(undef, fnw)
    iv = Array{ComplexF64}(undef, bnw)

    # Create arrays for storing functions
    for i in 1:fnw, j in 1:nx, k in 1:ny, l in 1:nz
        iw[i] = valueim(mesh.IR_basis_set.smpl_wn_f.sampling_points[i], beta) 
        Gw[i, j, k, l] = G(iw[i])
        Sigmaw[i, j, k, l] = Sigma(iw[i])
    end
    for i in 1:bnw, j in 1:nx, k in 1:ny, l in 1:nz
        iv[i] = valueim(mesh.IR_basis_set.smpl_wn_b.sampling_points[i], beta) 
        Vw[i, j, k, l] = V(iv[i])
    end

    # Evaluate the convolution
    Gt = wn_to_tau(mesh, Fermionic(), Gw)
    Vt = wn_to_tau(mesh, Bosonic(), Vw)
    Sigmat = Vt .* Gt 
    Sigmatw = tau_to_wn(mesh, Fermionic(), Sigmat)

    # Plot the results
    p = plot(real.(Sigmatw[:,1, 1, 1]), label="Calculated")
    plot!(p, imag.(Sigmatw[:, 1, 1, 1]), label="Calculated")
    plot!(p, real.(Sigmaw[:, 1, 1, 1]), label="Exact")
    plot!(p, imag.(Sigmaw[:, 1, 1, 1]), label="Exact")
    display(p)
    readline()
end

function test_ir_multiply()
    IR_basis_set = FiniteTempBasisSet(beta, 10.0, 1e-10)
    mesh = IRMesh.Mesh(IR_basis_set)
    fnw, bnw = mesh.fnw, mesh.bnw

    # Define constants
    g = 1.0
    wq = 1.0
    e = 0.4

    # Define functions
    V(w) = -g^2 * 2 * wq / (-wq^2 + w^2)
    G(w) = 1 / (w - e)
    n(e) = 1 / (exp(beta * e) - 1)
    f(e) = 1 / (exp(beta * e) + 1)
    Sigma(w) = g^2 * ( (1 + n(wq) - f(e)) / (w - (e + wq)) + (n(wq) + f(e)) / (w - (e - wq)))

    # Allocate arrays for storing functions
    Vw = Array{ComplexF64}(undef, bnw)
    Gw = Array{ComplexF64}(undef, fnw)
    Sigmaw = Array{ComplexF64}(undef, fnw)
    iw = Array{ComplexF64}(undef, fnw)
    iv = Array{ComplexF64}(undef, bnw)

    # Create arrays for storing functions
    for i in 1:fnw 
        iw[i] = valueim(mesh.IR_basis_set.smpl_wn_f.sampling_points[i], beta) 
        Gw[i] = G(iw[i])
        Sigmaw[i] = Sigma(iw[i])
    end
    for i in 1:bnw 
        iv[i] = valueim(mesh.IR_basis_set.smpl_wn_b.sampling_points[i], beta) 
        Vw[i] = V(iv[i])
    end

    # Evaluate the convolution
    Gt = wn_to_tau(mesh, Fermionic(), Gw)
    Vt = wn_to_tau(mesh, Bosonic(), Vw)
    Sigmat = Vt .* Gt 
    Sigmatw = tau_to_wn(mesh, Fermionic(), Sigmat)

    # Plot the results
    p = plot(real.(Sigmatw), label="Calculated")
    plot!(p, imag.(Sigmatw), label="Calculated")
    plot!(p, real.(Sigmaw), label="Exact")
    plot!(p, imag.(Sigmaw), label="Exact")
    display(p)
    readline()
end

function test_ir_multiply2()
    IR_basis_set = FiniteTempBasisSet(beta, 10.0, 1e-10)
    mesh = IRMesh.Mesh(IR_basis_set, nx, ny, nz)
    fnw, bnw, fntau, bntau = mesh.fnw, mesh.bnw, mesh.fntau, mesh.bntau

    # Define constants
    g = 1.0
    wq = 1.0

    # Define functions
    X(w) = -g^2 * 2 * wq / (-wq^2 + w^2)
    G(w, ek) = 1 / (w - ek)
    f(e) = 1 / (exp(beta * e) + 1)

    # Allocate arrays for storing functions
    Xw = Array{ComplexF64}(undef, bnw, nx, ny, nz)
    Gw = Array{ComplexF64}(undef, fnw, nx, ny, nz)
    iw = Array{ComplexF64}(undef, fnw)
    iv = Array{ComplexF64}(undef, bnw)

    # Create arrays for storing functions
    for i in 1:fnw, j in 1:nx, k in 1:ny, l in 1:nz
        iw[i] = valueim(mesh.IR_basis_set.smpl_wn_f.sampling_points[i], beta) 
        Gw[i, j, k, l] = G(iw[i], epsilon(get_kvec(j, k, l)))
    end
    for i in 1:bnw 
        iv[i] = valueim(mesh.IR_basis_set.smpl_wn_b.sampling_points[i], beta) 
        Xw[i, :, :, :] .= X(iv[i])
    end

    # Evaluate the convolution
    Gt = kw_to_rtau(Gw, 'F', mesh)
    Xt = Array{ComplexF64}(undef, bntau, nx, ny, nz)
    for i in 1:bntau, j in 1:nx, k in 1:ny, l in 1:nz
        Xt[i, j, k, l] = Gt[i, j, k, l] * Gt[bntau+1-i, j, k, l]
    end
    Xtw = rtau_to_kw(Xt, 'B', mesh)

    calc = Array{ComplexF64}(undef, nx)
    for i in 1:nx
        ind = trunc(Int, (bnw+1)/2)
        calc[i] = Xtw[ind, i, i, i]
    end
    # Plot the results
    p = plot(real.(calc), label="Calculated")
    display(p)
    readline()
end

function test_multiply()
    pts = 1000
    f = Array{ComplexF64}(undef, pts)
    exact_r = Array{ComplexF64}(undef, pts)
    IR_basis_set = FiniteTempBasisSet(beta, 10, 1e-10)
    mesh = IRMesh.Mesh(IR_basis_set, pts)
    for i in 1:pts
        x = 2 * pi * i / pts
        f[i] = sin(x)
        exact_r[i] = -cos(x) / 2 - sin(x) / (8 * pi)
    end
    temp = fft(f)
    r_r = temp .* temp
    calc_r = ifft(r_r) / pts

    p = plot(real.(calc_r), label="Calculated")
    plot!(p, real.(exact_r), label="Exact")
    display(p)
    readline()
end

function test_multiply_flat()
    pts = 1000
    f = Array{ComplexF64}(undef, pts)
    g = Array{ComplexF64}(undef, pts)
    exact_r = Array{ComplexF64}(undef, pts)
    IR_basis_set = FiniteTempBasisSet(beta, 10, 1e-10)
    mesh = IRMesh.Mesh(IR_basis_set, pts)
    sum = 0
    for i in 1:pts
        x = 1 * pi * i / pts
        f[i] = 1.0
        g[i] = sin(x)
        sum += f[i] * g[i]
    end
    exact_r .= sum / pts
    f_r = fft(f)
    g_r = fft(g)
    r_r = f_r .* g_r
    calc_r = ifft(r_r) / pts

    p = plot(real.(calc_r), label="Calculated")
    plot!(p, real.(exact_r), label="Exact")
    display(p)
    readline()
end

function test_multiply2()
    pts = 1000
    IR_basis_set = FiniteTempBasisSet(beta, 10, 1e-10)
    mesh = IRMesh.Mesh(IR_basis_set, pts)
    fnw, bnw = mesh.fnw, mesh.bnw
    f = Array{ComplexF64}(undef, fnw, pts)
    exact_r = Array{ComplexF64}(undef, fnw, pts)
    for i in 1:fnw, j in 1:pts
        a = 1 / abs(j - pts / 2)
        x = 2 * pi * j / pts
        f[i, j] = a * sin(x)
        exact_r[i, j] = a * (-cos(x) / 2 - sin(x) / (8 * pi))
    end
    temp = kw_to_rtau(f, 'F', mesh)
    r_r = temp .* temp
    calc_r = rtau_to_kw(r_r, 'F', mesh)

    p = plot(real.(calc_r[1,:]), label="Calculated")
    #plot!(p, real.(exact_r[1,:]), label="Exact")
    display(p)
    readline()
end

function test_transform()
    pts = 1000
    f = Array{ComplexF64}(undef, pts)
    g = Array{ComplexF64}(undef, pts)
    IR_basis_set = FiniteTempBasisSet(beta, 10, 1e-10)
    mesh = IRMesh.Mesh(IR_basis_set, pts)
    for i in 1:pts
        x = 2 * pi * i / pts
        f[i] = sin(x)
        g[i] = sin(x)
    end
    temp = fft(f)
    new_f = ifft(temp)

    p = plot(real.(new_f), label="Calculated")
    plot!(p, real.(g), label="Exact")
    display(p)
    readline()
end

function test_max_bcs_no_w()
    f = Array{ComplexF64}(undef, nx, ny)
    g = Array{ComplexF64}(undef, nx, ny)

    lambda = 1.0
    wD = 0.2
    mu = 1.0
    phi = 1.0
    exact = lambda * asinh(wD / phi) / (2 * pi)

    sphere(kx, ky) = kx^2 + ky^2
    func(e, z) = 1 / (e^2 + phi^2)^(0.5)
    almost_constant(e, z) = abs(e) > wD ? 0.0 : lambda
    for i in 1:nx, j in 1:ny
        kx = 2 * pi * (i / nx - 0.5)
        ky = 2 * pi * (j / ny - 0.5)
        e = sphere(kx, ky) - mu
        g[i, j] = func(e, 0)
        f[i, j] = almost_constant(e, 0)
    end

    ft = fft(f)
    gt = fft(g)
    resultt = ft .* gt
    result = ifft(resultt) / (nx * ny)

    println(maximum(real.(result)), " ", maximum(real.(exact)))
end


function test_max_bcs()
    IR_basis_set = FiniteTempBasisSet(beta, 10, 1e-10)
    mesh = IRMesh.Mesh(IR_basis_set, nx, ny)
    fnw, bnw = mesh.fnw, mesh.bnw
    f = Array{ComplexF64}(undef, bnw, nx, ny)
    g = Array{ComplexF64}(undef, fnw, nx, ny)
    exact = Array{ComplexF64}(undef, fnw)

    lambda = 1.0
    wD = 0.1
    mu = 1.0
    phi = 1.0

    sphere(kx, ky) = kx^2 + ky^2
    func(e, z) = 1 / (e^2 + phi^2 - z^2)
    almost_constant(e, z) = abs(e) > wD ? (lambda / (1 - 0.01 * z^2)) / 100 : lambda / (1 - 0.01 * z^2)
    #almost_constant(z) = lambda / (1 - 0.01 * z^2)

    for i in 1:fnw, j in 1:nx, k in 1:ny
        iw = valueim(mesh.IR_basis_set.smpl_wn_f.sampling_points[i], beta) 
        kx = 2 * pi * (j / nx - 0.5)
        ky = 2 * pi * (k / ny - 0.5)
        e = sphere(kx, ky) - mu
        g[i, j, k] = func(e, iw)
        #exact[i] = lambda / (2 * (a^2 + phi^2)^(0.5)) * tanh(beta * (a^2 + phi^2)^(0.5) / 2)
        exact[i] = lambda * asinh(wD / phi)
    end
    for i in 1:bnw, j in 1:nx, k in 1:ny
        iv = valueim(mesh.IR_basis_set.smpl_wn_b.sampling_points[i], beta)
        kx = 2 * pi * (j / nx - 0.5)
        ky = 2 * pi * (k / ny - 0.5)
        e = sphere(kx, ky) - mu
        f[i, j, k] = pi * almost_constant(e, iv)
    end

    ft = kw_to_rtau(f, 'B', mesh)
    gt = kw_to_rtau(g, 'F', mesh)
    resultt = ft .* gt
    result = rtau_to_kw(resultt, 'F', mesh)

    println(maximum(real.(result)), " ", maximum(real.(exact)))

    result_plot = sum(result, dims=(2, 3))[:, 1, 1] / (nx * ny)
    println("size: ", size(result_plot))
    p = plot(real.(result_plot), label="Calculated(Real)")
    plot!(imag.(result_plot), label="Calculated(Imag)")
    plot!(p, real.(exact), label="Exact")
    display(p)
    readline()
end

end # module

#using .Eliashberg
#Eliashberg.test_transform()
#Eliashberg.test_multiply_flat()
#Eliashberg.test_ir_multiply2()
#Eliashberg.gpu_sum()
#Eliashberg.eliashberg_sparse_ir()
#Eliashberg.test_max_bcs()
Eliashberg.test_max_bcs_no_w()
