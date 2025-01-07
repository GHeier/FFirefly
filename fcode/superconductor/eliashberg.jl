module Eliashberg

include("../objects/fields.jl")
using .Fields

using LinearAlgebra, Printf, PyCall
fcode = pyimport("fcode")
cfg = fcode.config

nx, ny, nz = cfg.k_mesh
println("nx: ", nx, " ny: ", ny, " nz: ", nz)
nw = cfg.w_pts
nw = 1
U = cfg.onsite_U

const beta = 1.0
const pi = π

function to_IBZ(k)
    q = abs.(k)
    q .= ifelse.(q .> π, -(q .- 2π), q)
    q .+= 1e-4
    return q
end

function V(k1, k2)
    qm = k1 - k2
    qm = to_IBZ(qm)
    Xm = chi(qm[1], qm[2], qm[3])
    Vm = U^2 * Xm / (1 - U*Xm) + U^3 * Xm^2 / (1 - U^2 * Xm^2)
    return Vm
end


function get_denominator(w1, phi_el, Z_el, chi_el, eps)
    return abs2(w1 * Z_el) + abs2(eps + chi_el) + abs2(phi_el)
end

function k_integral(k, w, w1, phi, Z, chi)
    phi_int = Z_int = chi_int = 0.0 + 0.0im
    n = round(Int, (imag(w) / (pi / beta) - 1.0) / 2.0)

    println("k: ", k)
    for i in 1:nx, j in 1:ny, l in 1:nz
        k1 = cfg.brillouin_zone * [i / nx, j / ny, l / nz]
        println("k1: ", k1)
        phi_el, Z_el, chi_el = phi[i, j, l, n], Z[i, j, l, n], chi[i, j, l, n]
        println("Phi: ", phi_el, " Z: ", Z_el, " Chi: ", chi_el)
        V_val = (V(k, k1) + V(k, -k1)) / 2.0
        denom = get_denominator(w1, phi_el, Z_el, chi_el, sum(abs2, k1))
        phi_int += V_val * phi_el / denom
        Z_int += V_val * w1 / w * Z_el / denom
        chi_int += V_val * (sum(abs2, k1) + chi_el) / denom
    end

    return phi_int / (2pi)^3, Z_int / (2pi)^3, chi_int / (2pi)^3
end

function eliashberg_sum(phi, Z, chi)
    new_phi = zeros(Complex{Float64}, size(phi))
    new_Z = zeros(Complex{Float64}, size(Z))
    new_chi = zeros(Complex{Float64}, size(chi))

    #Threads.@threads for a in 0:nx, b in 0:ny, c in 0:nz
    for a in 1:nx, b in 1:ny, c in 1:nz
        println("a: ", a, " b: ", b, " c: ", c)
        for i in 1:nw, j in 1:nw
            w = Complex(0.0, (2i+1)*pi / beta)
            w1 = Complex(0.0, (2j+1)*pi / beta + 0.0001)
            println("w: ", w, " w0: ", w1)
            k = cfg.brillouin_zone * [a / nx, b / ny, c / nz]
            phi_int, Z_int, chi_int = k_integral(k, w, w1, phi, Z, chi)
            new_phi[a, b, c, i] += phi_int / beta
            new_Z[a, b, c, i] += (1 + Z_int / beta)
            new_chi[a, b, c, i] -= chi_int / beta
        end
    end
    return new_phi, new_Z, new_chi
end

function evaluate_eliashberg()
    phi = Z = chi = ones(Complex{Float64}, nx, ny, nz, nw)
    iterations = 10

    input_data_file = "/home/g/Research/fcode/chi_mesh_dynamic.dat"
    global chi = Fields.create_interpolation(input_data_file, dims=cfg.dimension, field_type=:scalar)
    #fcode.load_global_chi(input_data_file)

    println("Starting Eliashberg calculation")
    for i in 0:iterations
        phi, Z, chi = eliashberg_sum(phi, Z, chi)
        println("(Sample) Phi: ", phi[1,1,1,1], " Z: ", Z[1,1,1,1], " Chi: ", chi[1,1,1,1])
    end
end

end # module

using .Eliashberg
Eliashberg.evaluate_eliashberg()
