module Eliashberg


include("../objects/fields.jl")
using .Fields
include("../objects/mesh.jl")
using .IRMesh

using SparseIR
import SparseIR: Statistics, value, valueim
using Base.Threads, JLD
using LinearAlgebra, Printf, PyCall
fcode = pyimport("fcode")
cfg = fcode.config

nx, ny, nz = cfg.k_mesh
dim = cfg.dimension
if dim == 2
    nz = 1
end
nw = cfg.w_pts 
U = cfg.onsite_U
BZ = cfg.brillouin_zone

const beta = 1/cfg.Temperature
println("Beta: ", beta)
const pi = π

function to_IBZ(k)
    tolerance = 1e-10  # small tolerance to account for floating-point errors
    q = abs.(k)
    q .= ifelse.(abs.(q .- 2π) .< tolerance, 0.0, ifelse.(q .> π, -(q .- 2π), q))
    q = abs.(q)
    q .+= 1e-4
    q .= ifelse.(q .> 2π, q .- 1e-4, q)
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
    Xm = Susceptibility(qm[1], qm[2], qm[3], iw)
    Vm = U^2 * Xm / (1 - U*Xm) + U^3 * Xm^2 / (1 - U^2 * Xm^2)
    return Vm
end

function get_bandwidth()
    maxval = -1000
    minval = 1000
    for i in 1:200, j in 1:200, k in 1:200
        kvec = BZ * [i / 200, j / 200, k / 200]
        eps = epsilon(kvec)
        maxval = max(maxval, eps)
        minval = min(minval, eps)
    end
    return maxval - minval
end

function fill_V_arr!(V_arr)
    Vw_arr = Array{Complex{Float64}}(undef, bnw, nx, ny)
    for i in 1:nx, j in 1:ny, k in 1:nz, l in 1:bnw
        kvec = BZ * [i / nx, j / ny, k / nz]
        w1 = Complex(0.0, pi*(2*l) / beta)
        Vw_arr[l, i, j] = V(kvec, [0,0,0], w1)
    end
    temp = k_to_r(mesh, Vw_arr)
    V_arr .= wn_to_tau(mesh, Bosonic(), temp)
    for i in 1:nx, j in 1:ny, k in 1:nz, l in 1:bntau
        V_arr[l, i, j] = V_arr[bntau - l + 1, i, j]
    end
    println("Filled V_arr")
end

function initialize_phi_Z!(phi_arr, Z_arr)
    for i in 1:inw, j in 1:nx, k in 1:ny
        w = pi*(2*(i-1) + 1) / beta
        kvec = BZ * [j / nx, k / ny, 0 / nz]
        dwave = (cos(kvec[1]) - cos(kvec[2])) / 2
        phi_arr[i, j, k] = 1 / (w + 1) * dwave
        Z_arr[i, j, k] = 1.0
    end
    println("Initialized phi and Z")
end

function condense_to_F!(phi, Z, F_arr, Fz_arr)
    for i in 1:nx, j in 1:ny, k in 1:nz, l in 1:inw
        w = Complex(0.0, pi*(2*(l-1) + 1) / beta)
        kvec = BZ * [i / nx, j / ny, k / nz]
        denom = get_denominator(w, phi[l, i, j], Z[l, i, j], epsilon(kvec))
        F_arr[l, i, j] = phi[l, i, j] / denom
        Fz_arr[l, i, j] = Z[l, i, j] * w / denom
    end
end

function update_phi_Z!(F_arr, Fz_arr, V_arr, phi_arr, Z_arr)
    temp = k_to_r(mesh, F_arr)
    F_arr .= wn_to_tau(mesh, Fermionic(), temp)
    temp = k_to_r(mesh, Fz_arr)
    Fz_arr .= wn_to_tau(mesh, Fermionic(), temp)

    phit = V_arr .* F_arr
    Zt = V_arr .* Fz_arr

    temp = r_to_k(mesh, phit)
    phi_arr .= tau_to_wn(mesh, Fermionic(), temp)
    temp = r_to_k(mesh, Z_arr)
    Z_arr .= tau_to_wn(mesh, Fermionic(), temp)

    for i in 1:inw
        w = Complex(0.0, pi*(2*(i-1) + 1) / beta)
        Z_arr[i, :, :] .= 1 .- 1.0 / (w*beta) .* Z_arr[i, :, :]
    end
    phi_arr .= phi_arr .* (-1.0 / beta)
end

function newfunc()
    println("Beginning Eliashberg")
    wmax = get_bandwidth()
    println("Bandwidth: ", wmax)
    IR_tol = 1e-10
    scf_tol = 1e-4
    IR_basis_set = FiniteTempBasisSet(beta, Float64(wmax), IR_tol)
    global mesh = Mesh(nx, ny, IR_basis_set)
    global inw = mesh.fnw
    global bnw = mesh.bnw
    global bntau = mesh.bntau
    println("IRMesh created")
    input_data_file = "/home/g/Research/fcode/chi_mesh_dynamic.dat"
    global Susceptibility = Fields.create_interpolation(input_data_file, dims=4, field_type=:scalar, value_type=ComplexF64)

    phi_arr = Array{ComplexF64}(undef, inw, nx, ny)
    Z_arr = Array{ComplexF64}(undef, inw, nx, ny)
    initialize_phi_Z!(phi_arr, Z_arr)

    println("maxphi: ", maximum(abs.(phi_arr)))
    println("maxZ: ", maximum(abs.(Z_arr)))

    V_arr = Array{ComplexF64}(undef, mesh.bntau, mesh.nk1, mesh.nk2)
    fill_V_arr!(V_arr)
    println("maxV: ", V_arr[argmax(abs.(V_arr))])
    println("maxV: ", maximum(abs.(V_arr)))

    F_arr = Array{ComplexF64}(undef, inw, nx, ny)
    Fz_arr = Array{ComplexF64}(undef, inw, nx, ny)
    condense_to_F!(phi_arr, Z_arr, F_arr, Fz_arr)
    println("maxF: ", F_arr[argmax(abs.(F_arr))])
    println("maxFz: ", Fz_arr[argmax(abs.(Fz_arr))])


    max_phi = max_Z = 0.0
    prev_phi = prev_Z = 1.0
    println("Starting Convergence Loop")
    for i in 1:100
        update_phi_Z!(F_arr, Fz_arr, V_arr, phi_arr, Z_arr)

        max_phi = maximum(abs.(phi_arr))
        max_Z = maximum(abs.(Z_arr))

        println("\n $i) Max phi: ", max_phi)
        println(" $i) Max Z: ", max_Z)
        error = abs(max_phi / prev_phi) - abs(max_Z / prev_Z)
        println(" $i) Error: ", error)
        if abs(error) < scf_tol
            println("Converged after ", i, " iterations")
            break
        end
        prev_phi, prev_Z = max_phi, max_Z
        condense_to_F!(phi_arr, Z_arr, F_arr, Fz_arr)
    end
end 

function get_denominator(w1, phi_el, Z_el, eps)
    return abs2(w1 * Z_el) + abs2(eps) + abs2(phi_el)
end


function k_integral(k, w, w1, phi, Z)
    phi_int = Z_int = 0.0 + 0.0im
    n = trunc(Int, (imag(w) * beta / pi - 1) / 2 + 1)
    temp = 0

    for i in 1:nx 
        for j in 1:ny 
            for l in 1:nz
                k1 = BZ * [i / nx, j / ny, l / nz]
                phi_el, Z_el= phi[i, j, l, n], Z[i, j, l, n]
                V_val = (V(k, k1, w) + V(k, -k1, w)) / 2.0
                denom = get_denominator(w1, phi_el, Z_el, epsilon(k1))
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

end # module

using .Eliashberg
println(Threads.nthreads(), " threads available")
#Eliashberg.evaluate_eliashberg()
Eliashberg.newfunc()
