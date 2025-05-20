module Linearized_Eliashberg
println("Started Julia")
using Firefly

include("../objects/mesh.jl")
using .IRMesh

using SparseIR
import SparseIR: Statistics, value, valueim, MatsubaraSampling64F, TauSampling64
using Roots, LinearAlgebra, Printf, DelimitedFiles
using FFTW
using Base.Threads, MPI

cfg = Firefly.Config

const prefix = cfg.prefix 
const outdir = cfg.outdir 

nx, ny, nz = cfg.k_mesh
const dim = cfg.dimension
if dim == 2 nz = 1 end
const nk = nx * ny * nz
const nw = cfg.w_pts
const nbnd = cfg.nbnd

const U = cfg.onsite_U
const BZ = cfg.brillouin_zone
const mu = cfg.fermi_energy
const wc = cfg.bcs_cutoff_frequency

const beta = 1 / cfg.Temperature

const bcs_debug = false 

if bcs_debug && nw % 2 != 0
    println("Remember to make nw even")
    exit()
end

function get_kvec(i, j, k)
    temp = BZ * [i / nx - 0.5, j / ny - 0.5, k / nz - 0.5]
    return temp[1:dim]
end

function square_well_approximate(w)
    delta = 0.01
    return - 0.5 * (tanh((w + wc)/delta) - tanh((w - wc)/delta))
end

function zero_out_beyond_wc!(A, iw)
    if size(A, 1) != length(iw)
        println("Zero_Out: Incorrect array mapping")
    end
    for i in 1:length(iw)
        if abs(imag(iw[i])) > wc 
            A[i, :, :, :] .= 0
        end
    end
end

function fill_V_rt(vertex, iv, mesh)
    bnw = length(iv)
    V = Array{ComplexF32}(undef, bnw, nx, ny, nz)
    for i in 1:nx, j in 1:ny, k in 1:nz, l in 1:bnw
        q = get_kvec(i, j, k)
        w = imag(iv[l])
        V[l, i, j, k] = 0.5 * (vertex(q, w) + vertex(-q, w))
    end
    return kw_to_rtau(V, 'B', mesh)
end


function get_iw_iv(mesh)
    fnw, bnw, fntau, bntau = mesh.fnw, mesh.bnw, mesh.fntau, mesh.bntau
    println("fnw: ", fnw, " bnw: ", bnw, " fntau: ", fntau, " bntau: ", bntau)
    println("Created IRMesh")
    iw = Array{ComplexF32}(undef, fnw)
    iv = Array{ComplexF32}(undef, bnw)
    for i in 1:fnw
        iw[i] = valueim(mesh.IR_basis_set.smpl_wn_f.sampling_points[i], beta)
    end
    for i in 1:bnw
        iv[i] = valueim(mesh.IR_basis_set.smpl_wn_b.sampling_points[i], beta)
    end
    println("Fermionic frequency from ", minimum(imag.(iw)), " to ", maximum(imag.(iw)))
    println("Bosonic frequency from ", minimum(imag.(iv)), " to ", maximum(imag.(iv)))
    return iw, iv 
end


function create_energy_mesh(band)
    e_arr = zeros(Float32, nx, ny, nz)
    for i in 1:nx, j in 1:ny, k in 1:nz
        e_arr[i, j, k] = band(1, get_kvec(i, j, k)) - mu
    end
    e_4d = reshape(e_arr, 1, nx, ny, nz)
    return e_4d
end


function Z_convergence!(Z, V_rt, iw, e, mesh)
    Z_err = 1.0
    Z_max_i = 0.0
    Z_max_f = 0.0
    while Z_err > 1e-4
        F = imag.(iw) .* Z ./ ( (imag.(iw) .* Z) .^2 .+ e.^2)
        F_rt = kw_to_rtau(F, 'F', mesh)

        temp = V_rt .* F_rt
        result = rtau_to_kw(temp, 'F', mesh)

        Z .= 1 .- result ./ imag.(iw)

        Z_max_f = maximum(abs.(real.(Z)))
        Z_err = abs(Z_max_f - Z_max_i)
        Z_max_i = Z_max_f
        println("Z_max: ", Z_max_f)
        println("Z_min: ", minimum(abs.(real.(Z))))
    end
end


function linearized_eliashberg_loop(phi, Z, iw, V_rt, e, mesh = 0)
    eig, prev_eig, eig_err = 0, 0, 1
    if bcs_debug
        N = Firefly.Field_R(outdir * prefix * "_DOS.dat")
        initial_expected = log(1.134 * wc * beta) * N(mu)
        println("Initial Expected eig: $initial_expected")
    end
    while eig_err > 1e-5
        if !bcs_debug
            eig, phi = linearized_eliashberg(phi, Z, iw, V_rt, e, mesh)
        else
            eig, phi = linearized_eliashberg_bcs(phi, Z, iw, V_rt, e)
        end
        eig_err = abs(eig - prev_eig)
        prev_eig = eig
        println("Eig: ", round(eig, digits=6), " Error: ", round(eig_err, digits=6), "   \r")
        return eig, phi
    end
    println()
    if bcs_debug
        final_expected = log(1.134 * wc * beta)
        println("Final Expected eig: $final_expected")
    end
    return eig, phi
end


function linearized_eliashberg(phi, Z, iw, V_rt, e, mesh)
    F = -phi ./ ((Z .* imag.(iw)).^2 .+ e.^2)
    F_rt = kw_to_rtau(F, 'F', mesh)
    temp = V_rt .* F_rt
    result = rtau_to_kw(temp, 'F', mesh)
    result = phi .* result # Resultant eigenvector/value
    mag = phi .* phi # Normalized gap

    eig = sum(result) / (nk * beta)
    norm = sum(mag) / (nk * beta)
    return eig / norm, result ./ norm
end


function linearized_eliashberg_bcs(phi, Z, iw, V_rt, e)
    zero_out_beyond_wc!(phi, iw)

    F = -phi ./ ((Z .* imag.(iw)).^2 .+ e.^2)
    #F .= -1.0
    F_rt = fft(F)
    temp = V_rt .* F_rt
    result = fftshift(ifft(temp)) / (beta * nk)
    zero_out_beyond_wc!(result, iw)

    result = phi .* result # Resultant eigenvector/value
    mag = phi .* phi # Normalized gap

    eig = sum(result) / (nk * beta)
    norm = sum(mag) / (nk * beta)
    return eig / norm, result ./ norm
end


function linear_index(i, j, k, l, nx, ny, nz, fnw)
    return l + fnw * (k - 1 + nz * (j - 1 + ny * (i - 1)))
end

function fill_matrix(F, iw, vertex)
    fnw = length(iw)
    len = nk * fnw
    println("Allocating $(round(len^2 * 16 / 1024^2, digits=2)) MB")
    matrix = Array{ComplexF32}(undef, len, len)

    for i in 1:nx, j in 1:ny, k in 1:nz, l in 1:fnw
        kvec1 = get_kvec(i, j, k)
        w1 = iw[l]
        ind1 = linear_index(i, j, k, l, nx, ny, nz, fnw)
        for a in 1:nx, b in 1:ny, c in 1:nz, d in 1:fnw
            kvec2 = get_kvec(a, b, c)
            w2 = iw[d]
            ind2 = linear_index(a, b, c, d, nx, ny, nz, fnw)
            matrix[ind1, ind2] = vertex(kvec1 - kvec2, w1 - w2) * F[d, a, b, c]
        end
    end
    return matrix
end


function linearized_eliashberg_diagonalize(Z, iw, vertex, e_4d, mesh)
    F = -1.0 ./ ((Z .* iw).^2 .+ e_4d.^2)
    matrix = fill_matrix(F, iw, vertex) 
    E = eigen(matrix)
    V = E.vectors
    println("First 5", E[:5])
    println("Last 5", E[end-5:end])
    return E, V
end


function save_gap(phi, iw)
    println("Saving Gap")
    nw = length(iw)
    filename = outdir * prefix * "_gap.dat"
    open(filename, "w") do io
        if dim == 3
            println(io, "kx\tky\tkz\tw\tRe(f)\tIm(f)")
        else 
            println(io, "kx\tky\tw\tRe(f)\tIm(f)")
        end


        for i in 1:nx, j in 1:ny, k in 1:nz, l in 1:nw
            kvec = get_kvec(i, j, k)
            w = iw[l]
            f = phi[l, i, j, k]
            if dim == 3
                kx, ky, kz = kvec
                println(io, "$(kx)\t$(ky)\t$(kz)\t$(w)\t$(real(f))\t$(imag(f))")
            else
                kx, ky = kvec
                println(io, "$(kx)\t$(ky)\t$(w)\t$(real(f))\t$(imag(f))")
            end
        end
    end
    println("Gap saved to ", filename)
end


function eigenvalue_computation()
    if !bcs_debug
        println("Creating Mesh")
        mesh = IR_Mesh()
        iw, iv = get_iw_iv(mesh)
        fnw, bnw = length(iw), length(iv)
    else 
        mesh = 0
        iw, iv = get_iw_iv_bcs()
        fnw, bnw = nw, nw
        println("BCS DEBUG: Found iw/iv")
    end


    if !bcs_debug
        println("Getting Vertex")
        vertex = Firefly.Vertex()
        println("Fourier Transforming Vertex")
        V_rt = fill_V_rt(vertex, iv, mesh)
    else
        V = -1.0 .* ones(Complex{Float32}, bnw, nx, ny, nz)
        V_rt = fft(V)
        println("BCS DEBUG: FFT'd Vertex")
    end

    println("Getting Bands")
    band = Firefly.Bands()
    println("Creating energy mesh")
    e_4d = create_energy_mesh(band)

    println("Calculating Z")
    Z = ones(Complex{Float32}, fnw, nx, ny, nz)
    if !bcs_debug
        Z_convergence!(Z, V_rt, iw, e_4d, mesh)
    end

    #println("Diagonalization to find Eigenvalues and symmetries")
    #eigs, phis = linearized_eliashberg_diagonalize(Z, iw, V_rt, e_4d, mesh)

    println("Power Iteration to find Eigenvalue and symmetry")
    #phi = rand(Float32, fnw, nx, ny, nz) .+ im * rand(Float32, fnw, nx, ny, nz)
    phi = ones(fnw, nx, ny, nz)
    eig, phi = linearized_eliashberg_loop(phi, Z, iw, V_rt, e_4d, mesh)
    if real(eig) < 0 
        V = rtau_to_kw(V_rt, 'F', mesh) .- eig
        V_rt = kw_to_rtau(V, 'F', mesh)
        eig, phi = linearized_eliashberg_loop(phi, Z, iw, V_rt, e_4d, mesh)
    end
    @printf("Final Eig: %.6f\n", real(eig))
    if !bcs_debug
        save_gap(phi, iw)
    end
end


function get_iw_iv_bcs()
    iw = Array{ComplexF32}(undef, nw)
    iv = Array{ComplexF32}(undef, nw)
    for i in 1:nw
        n = -nw / 2 + i
        iw[i] = pi / beta * (2 * n + 1) * im
    end
    for i in 1:nw
        n = -nw / 2 + i
        iv[i] = pi / beta * (2 * n) * im
    end
    println("Fermionic frequency from ", minimum(imag.(iw)), " to ", maximum(imag.(iw)))
    println("Bosonic frequency from ", minimum(imag.(iv)), " to ", maximum(imag.(iv)))
    return iw, iv 
end

function square_well_test()
    N = Firefly.Field_R(outdir * prefix * "_DOS.dat")
    T = 1 / beta
    nk = nx * ny * nz
    epts = 1000
    de = 8 / epts
    println("Max n: ", (wc / (pi * T) - 1) / 2)
    phi = ones(ComplexF64, nw)
    result_phi = zeros(ComplexF64, nw)
    for k in 1:nw
        n = -nw / 2 + k
        iw_n = pi * T * (2 * n + 1)
        if abs(iw_n) > wc
            phi[k] = 0
        end
    end

    for k in 1:nw
        sum = 0
        n = -nw / 2 + k
        iw_n = pi * T * (2 * n + 1) * im
        for i in 1:epts, j in 1:nw
            e = -4 + (i - 1) * de
            m = -nw / 2 + j
            iw_m = pi * T * (2 * m + 1) * im
            V = 1
            #V = square_well_approximate(imag(iw_n - iw_m))
            if abs(imag(iw_m)) > wc
                V = 0
            end
            sum += N(e) * V * phi[j] / (imag(iw_m)^2 + (e - mu)^2)
        end
        sum *= T * de
        result_phi[k] = sum
    end
    println("Max: ", maximum(real.(result_phi)), " Min: ", minimum(real.(result_phi)))
    writedlm("phi_nw.dat", real.(result_phi))
end

end # module
