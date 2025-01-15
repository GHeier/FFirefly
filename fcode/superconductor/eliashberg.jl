module Eliashberg

include("../objects/fields.jl")
using .Fields
include("../objects/mesh.jl")
using .Mesh

using Base.Threads, JLD
using LinearAlgebra, Printf, PyCall
fcode = pyimport("fcode")
cfg = fcode.config

nx, ny, nz = cfg.k_mesh
dim = cfg.dimension
println("dim: ", dim)
if dim == 2
    nz = 1
end
println("nx: ", nx, " ny: ", ny, " nz: ", nz)
nw = cfg.w_pts 
println("nw: ", nw)
U = cfg.onsite_U
BZ = cfg.brillouin_zone
println("typeof(BZ): ", typeof(BZ))

const beta = 1/cfg.Temperature
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
            @inbounds begin
                new_phi[a, b, c, i] -= phi_sum / beta
                new_Z[a, b, c, i] += (1 - Z_sum / beta)
                #println("phi_sum: ", phi_sum, " Z_sum: ", Z_sum)
            end
        end
    end
    return new_phi, new_Z
end

function evaluate_eliashberg()
    phi = Z = ones(Complex{Float64}, nx, ny, nz, nw)
    iterations = 100

    println(Threads.nthreads(), " threads available")
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
Eliashberg.evaluate_eliashberg()
