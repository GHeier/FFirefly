module Eliashberg

using LinearAlgebra, Printf, PyCall

fcode = pyimport("fcode")

nx, ny, nz = fcode.k_mesh
nw = fcode.w_pts

const beta = 1.0
const pi = π

function get_denominator(w1, phi_el, Z_el, chi_el, eps)
    return abs2(w1 * Z_el) + abs2(eps + chi_el) + abs2(phi_el)
end

function k_integral(k, w, w1, phi, Z, chi)
    phi_int = Z_int = chi_int = 0.0 + 0.0im
    n = round(Int, (imag(w) / (pi / beta) - 1.0) / 2.0)

    for i in 1:nx, j in 1:ny, l in 1:nz
        k1 = [-pi + i * 2pi / nx, -pi + j * 2pi / ny, -pi + l * 2pi / nz]
        phi_el, Z_el, chi_el = phi[i, j, l, n], Z[i, j, l, n], chi[i, j, l, n]
        V = (arr[i, j, l] + arr[i, j, l]) / 2
        denom = get_denominator(w1, phi_el, Z_el, chi_el, sum(abs2, k1))
        phi_int += V * phi_el / denom
        Z_int += V * w1 / w * Z_el / denom
        chi_int += V * (sum(abs2, k1) + chi_el) / denom
    end

    return phi_int / (2pi)^3, Z_int / (2pi)^3, chi_int / (2pi)^3
end

function eliashberg_sum(phi, Z, chi)
    new_phi = zeros(Complex{Float64}, size(phi))
    new_Z = zeros(Complex{Float64}, size(Z))
    new_chi = zeros(Complex{Float64}, size(chi))

    Threads.@threads for a in 1:nx, b in 1:ny, c in 1:nz
        for i in 1:nw
            w = Complex(0.0, (2i+1)*pi / beta)
            for j in 1:nw
                w1 = Complex(0.0, (2j+1)*pi / beta + 0.0001)
                phi_int, Z_int, chi_int = k_integral([a, b, c], w, w1, phi, Z, chi)
                new_phi[a, b, c, i] += phi_int / beta
                new_Z[a, b, c, i] += (1 + Z_int / beta)
                new_chi[a, b, c, i] -= chi_int / beta
            end
        end
    end
    return new_phi, new_Z, new_chi
end

function evaluate_eliashberg()
    phi = Z = chi = ones(Complex{Float64}, nx, ny, nz, nw)
    iterations = 10
    susceptibility = fcode.ComplexField("chi_static_mesh.dat")

    println("Starting Eliashberg calculation")
    for i in 1:iterations
        phi, Z, chi = eliashberg_sum(phi, Z, chi)
        println("(Sample) Phi: ", phi[1,1,1,1], " Z: ", Z[1,1,1,1], " Chi: ", chi[1,1,1,1])
    end
end

end # module

