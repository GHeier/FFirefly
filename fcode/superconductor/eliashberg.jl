module Eliashberg

t1 = time()
include("../objects/field.jl")
using .CondensedMatterField
include("../objects/mesh.jl")
using .IRMesh

using CUDA, FFTW
using Plots
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
if dim == 2
    nz = 1
end
nw = cfg.w_pts 
U = cfg.onsite_U
BZ = cfg.brillouin_zone
mu = cfg.fermi_energy
t4 = time()
println("Time to load config: ", t4 - t3)

const beta = 1/cfg.Temperature
const pi = π



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
    n = (beta * imag(w) / pi - 1) / 2
    lambda = 2.0
    v = beta * 1.5 / (2 * pi)
    return -2 * lambda * v^2 / (n^2 + v^2)
end

function fill_V_arr!(V_arr, iw)
    Vw_arr = Array{Complex{Float64}}(undef, bnw, nx, ny)
    for i in 1:nx, j in 1:ny, k in 1:nz, l in 1:bnw
        kvec = get_kvec(i, j, k)
        w1 = iw[l]
        #Vw_arr[l, i, j] = V(kvec, zeros(dim), w1)
        Vw_arr[l, i, j] = paper2_V(w1)
    end
    temp = k_to_r(Vw_arr)
    V_arr .= wn_to_tau(mesh, Bosonic(), temp)
    println("Filled V_arr")
end

function initialize_phi_Z_chi!(phi_arr, Z_arr, chi_arr, iw)
    for i in 1:fnw, j in 1:nx, k in 1:ny, l in 1:nz
        w = iw[i]
        kvec = get_kvec(j, k, l)
        dwave = (cos(kvec[1]) - cos(kvec[2])) / 2
        phi_arr[i, j, k] = 1 / abs(imag(w)^2 + 1)# * dwave
        Z_arr[i, j, k] = 1.0 
        chi_arr[i, j, k] = 0.0
    end
end

function get_kvec(i, j, k)
    temp = BZ * [i / nx - 0.5, j / ny - 0.5, k / nz - 0.5]
    return temp[1:dim]
end

function condense_to_F_and_G!(phi, Z, chi, F_arr, G_arr, iw)
    for i in 1:nx, j in 1:ny, k in 1:nz, l in 1:fnw
        w = iw[l]
        kvec = get_kvec(i, j, k)
        denom = get_denominator(w, phi[l, i, j], Z[l, i, j], chi[l, i, j] + epsilon(kvec) - mu)
        F_arr[l, i, j] = -phi[l, i, j] / denom
        G_arr[l, i, j] = -(w * Z[l, i, j] + epsilon(kvec) - mu + chi[l, i, j]) / denom
    end
end

function update!(F, G, V_arr, phi, Z, chi, iw, sigma)
    F_rt = kw_to_rtau(F, 'F', mesh)
    G_rt = kw_to_rtau(G, 'F', mesh)

    phit = V_arr .* F_rt / mesh.nk
    sigmat = -V_arr .* G_rt / mesh.nk

    phi .= rtau_to_kw(phit, 'F', mesh)
    sigma .= rtau_to_kw(sigmat, 'F', mesh)

    for i in 1:fnw, j in 1:nx, k in 1:ny
        w = iw[i]
        sp, sm = sigma[i, j, k], sigma[fnw - i + 1, j, k]
        Z[i, j, k] = 1.0 - 0.5 * (sp - sm) / w
        chi[i, j, k] = 0.5 * (sp + sm)
    end
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
    global Susceptibility = CondensedMatterField.CMF(input_data_file)
    println("Loaded Susceptibility")

    phi_arr = Array{ComplexF64}(undef, fnw, nx, ny)
    Z_arr = Array{ComplexF64}(undef, fnw, nx, ny)
    chi_arr = Array{ComplexF64}(undef, fnw, nx, ny)
    println("Initializing phi, Z, and chi")
    initialize_phi_Z_chi!(phi_arr, Z_arr, chi_arr, iw)

    println("Initializing V")
    sigma = Array{ComplexF64}(undef, fnw, nx, ny, nz)
    V_arr = Array{ComplexF64}(undef, bntau, nx, ny, nz)
    fill_V_arr!(V_arr, iv)

    println("Initializing F and G")
    F_arr = Array{ComplexF64}(undef, fnw, nx, ny, nz)
    G_arr = Array{ComplexF64}(undef, fnw, nx, ny, nz)
    condense_to_F_and_G!(phi_arr, Z_arr, chi_arr, F_arr, G_arr, iw)


    phierr = 0.0
    prev_phi_arr = copy(phi_arr)
    println("Starting Convergence Loop")
    iterations = 300
    for i in 1:iterations
        update!(F_arr, G_arr, V_arr, phi_arr, Z_arr, chi_arr, iw, sigma)

        phierr = round(sum(abs.(prev_phi_arr - phi_arr)) / (fnw * nx * ny),digits=8)
        print("Iteration $i: Error = $phierr           \r")

        if abs(phierr) < scf_tol || maximum(abs.(phi_arr)) < 1e-10
            println("Converged after ", i, " iterations              ")
            break
        end
        prev_phi_arr = copy(phi_arr)
        condense_to_F_and_G!(phi_arr, Z_arr, chi_arr, F_arr, G_arr, iw)
    end

    max_phi = maximum(abs.(phi_arr))
    max_Z = maximum(abs.(Z_arr))

    println("Error: ", phierr)
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
    println("Saving phi and Z")
    save_arr_as_meshgrid(points, reordered_phi, "phi.dat", true)
    save_arr_as_meshgrid(points, reordered_Z, "Z.dat", true)
end 

function get_denominator(w1, phi_el, Z_el, eps)
    return abs2(w1 * Z_el) + abs2(eps) + abs2(phi_el)
end


function k_integral(k, w, w1, phi, Z)
    phi_int = Z_int = 0.0 + 0.0im
    n = trunc(Int, (imag(w) * beta / pi - 1) / 2)
    temp = 0

    for i in 1:nx 
        for j in 1:ny 
            for l in 1:nz
                k1 = BZ * [i / nx, j / ny, l / nz]
                phi_el, Z_el= phi[i, j, l, n], Z[i, j, l, n]
                V_val = (V(k, k1, w) + V(k, -k1, w)) / 2.0
                denom = get_denominator(w1, phi_el, Z_el, epsilon(k1) - mu)
                temp = max(temp, denom)
                phi_int += V_val * phi_el / denom
                Z_int += V_val * w1 / w * Z_el / denom
            end
        end
    end

    #println(w.im, " ", w1.im, " ", Z_int.re, " ", temp)
    return phi_int / (nx * ny * nz), Z_int / (nx * ny * nz)
end

function eliashberg_sum(phi, Z)
    new_phi = zeros(Complex{Float64}, size(phi))
    new_Z = zeros(Complex{Float64}, size(Z))

    for a in 1:nx, b in 1:ny, c in 1:nz
        print("Iteration ", (a-1)*nx^(dim-1) + (b-1)*ny^(dim-2) + c-1, " out of ", nx*ny*nz, "      \r")
        flush(stdout)
        @threads for i in 1:nw 
            phi_sum = Z_sum = 0.0 + 0.0im
            for j in 1:nw
                w = Complex(0.0, pi*(2*(i-1) + 1) / beta)
                w1 = Complex(0.0, pi*(2*(j-1) + 1) / beta)
                k = BZ * [a / nx, b / ny, c / nz]
                phi_int, Z_int = k_integral(k, w, w1, phi, Z)
                phi_sum += phi_int
                Z_sum += Z_int
            end
            new_phi[a, b, c, i] -= phi_sum / beta
            new_Z[a, b, c, i] += (1 - Z_sum / beta)
            #println("phi_sum: ", phi_sum, " Z_sum: ", Z_sum)
        end
    end
    return new_phi, new_Z
end

function evaluate_eliashberg()
    phi = Z = ones(Complex{Float64}, nx, ny, nz, nw)
    iterations = 100

    input_data_file = "/home/g/Research/fcode/chi_mesh_dynamic.dat"
    global Susceptibility = Fields.create_interpolation(input_data_file, dims=4, field_type=:scalar, value_type=ComplexF64)
    value = Susceptibility(0.005, 0.005, 0.005, 0.005)
    println("Value: ", value)
    #fcode.load_global_chi(input_data_file)

    println("Starting Eliashberg calculation")
    prev_phi = prev_Z = 0.0
    for i in 0:iterations
        phi, Z = eliashberg_sum(phi, Z)
        max_phi = maximum(abs.(phi))
        max_Z = maximum(abs.(Z))
        println("\nMax phi: ", max_phi)
        println("Max Z: ", max_Z)
        error = abs(prev_phi - max_phi) + abs(prev_Z - max_Z)
        println("Error: ", maximum(error))
        if maximum(error) < 1e-4
            println("Converged after ", i, " iterations")
            break
        end
        prev_phi, prev_Z = max_phi, max_Z
    end
    save("phi.jld", "phi", phi)
    save("Z.jld", "Z", Z)
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
        npts[i] = -fnw/2 + i - 1
        iw[i] = Complex(0.0, (pi*(2*(npts[i] - 1) + 1) / beta))
        phi[i] = 1 / abs(iw[i] + 1)
        Z[i] = 1 / abs(iw[i] + 1)
    end

    iw_gpu = CuArray(iw)
    phi_gpu = CuArray(phi)
    Z_gpu = CuArray(Z)

    BZ_gpu = CuArray(BZ)

    for iters in 1:20
        print("Iteration $iters: ")
        new_phi = CuArray(zeros(Complex{Float64}, fnw))
        new_Z = CuArray(ones(Complex{Float64}, fnw))

        #@device_code_warntype kernel_summation!(
        #    new_phi, new_Z, phi_gpu, Z_gpu, iw_gpu, nx, ny, beta, BZ_gpu, mu
        #)
        # GPU Kernel for summation
        @cuda threads=256 blocks=cld(fnw, 256) kernel_summation!(
             new_phi, new_Z, phi_gpu, Z_gpu, iw_gpu, nx, ny, beta, BZ_gpu, mu
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
    beta::Float64,
    BZ::CuDeviceMatrix{Float64},
    mu::Float64
)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    if i > length(new_phi)
        return
    end

    temp_phi = ComplexF64(0.0, 0.0)
    temp_Z = ComplexF64(0.0, 0.0)

    V(x) = 1.0 / (imag(x)^2 + 0.01)
    e(kx, ky, kz) = -2 * (cos(kx) + cos(ky))
    get_denom(iw, phi, Z, eps) = abs2(imag(iw)*Z) + abs2(eps) + abs2(phi)

    for j in 1:length(phi)
        for k in 1:nx, l in 1:ny
            kx = BZ[1, 1] * (k / nx) + BZ[1, 2] * (l / ny) + BZ[1, 3] * 0.0
            ky = BZ[2, 1] * (k / nx) + BZ[2, 2] * (l / ny) + BZ[2, 3] * 0.0
            kz = BZ[3, 1] * (k / nx) + BZ[3, 2] * (l / ny) + BZ[3, 3] * 0.0
            denom = get_denom(iw[j], phi[j], Z[j], e(kx, ky, kz) - mu)
            w = iw[i] - iw[j]

            V_val = ComplexF64(-1.0, 0.0)
            V_val = V(w)
            # Accumulate results
            temp_phi -= (1 / beta) * V_val * phi[j] / denom
            temp_Z -= (1 / beta) * V_val * (iw[j] / iw[i]) * Z[j] / denom
        end
    end

    new_phi[i] = temp_phi
    new_Z[i] = temp_Z
    return nothing
end

end

using .Eliashberg
Eliashberg.eliashberg_sparse_ir()
