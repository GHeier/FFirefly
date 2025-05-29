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

function zero_out_beyond_wc_e!(A, e)
    nw, nx, ny, nz = size(A)
    for j in 1:nx, k in 1:ny, l in 1:nz
        if abs(e[1, j, k, l]) > wc 
            A[:, j, k, l] .= 0
        end
    end
end

function zero_out_beyond_wc_iw!(A, iw)
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
        V[l, i, j, k] = vertex(q, w)
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


function create_energy_mesh(band, iw, Sigma)
    nw = length(iw)
    e_arr = zeros(ComplexF32, nw, nx, ny, nz)
    for i in 1:nx, j in 1:ny, k in 1:nz, l in 1:nw
        kvec = get_kvec(i, j, k)
        e_arr[l, i, j, k] = band(1, kvec) - mu + Sigma(kvec, imag(iw[l]))
    end
    return e_arr
end


function linearized_eliashberg_loop(phi, iw, V_rt, e, mesh = 0)
    eig, prev_eig, eig_err = 0, 0, 1
    if bcs_debug
        N = Firefly.Field_R(outdir * prefix * "_DOS.dat")
        initial_expected = log(1.134 * wc * beta) * N(mu)
        println("Initial Expected eig: $initial_expected")
    end
    while eig_err > 1e-5
        if !bcs_debug
            eig, phi = linearized_eliashberg(phi, iw, V_rt, e, mesh)
        else
            eig, phi = linearized_eliashberg_bcs(phi, iw, V_rt, e)
        end
        eig_err = abs(eig - prev_eig)
        prev_eig = eig
        println("Eig: ", round(eig, digits=6), " Error: ", round(eig_err, digits=6), "   \n")
    end
    println()
    if bcs_debug
        final_expected = log(1.134 * wc * beta)
        println("Final Expected eig: $final_expected")
    end
    return eig, phi
end


function linearized_eliashberg(phi, iw, V_rt, e, mesh)
    F = -phi .* ( 1 ./ (iw .- e)) .* (1 ./ (-iw .- e))
    #F = -phi ./ ((Z .* imag.(iw)).^2 .+ e.^2)
    F_rt = kw_to_rtau(F, 'F', mesh)
    temp = V_rt .* F_rt
    newphi = rtau_to_kw(temp, 'F', mesh)

    mag = newphi .* newphi # Normalized gap
    norm = sum(mag)^(0.5)

    result = phi .* newphi # Resultant eigenvector/value
    eig = sum(result)
    return eig, newphi ./ norm
end


function linearized_eliashberg_bcs(phi, iw, V_rt, e)
    zero_out_beyond_wc_e!(phi, e)
    #zero_out_beyond_wc_iw!(phi, iw)

    F = -phi .* ( 1 ./ (iw .- e)) .* (1 ./ (-iw .- e))
    #F = -phi ./ ((Z .* imag.(iw)).^2 .+ e.^2)
    #F .= -1.0
    F_rt = fft(F)
    temp = V_rt .* F_rt
    result = fftshift(ifft(temp)) / (beta * nk)
    zero_out_beyond_wc_e!(phi, e)
    #zero_out_beyond_wc_iw!(result, iw)

    result = phi .* result # Resultant eigenvector/value
    mag = phi .* phi # Normalized gap

    eig = sum(result) / (nk * beta)
    norm = sum(mag) / (nk * beta)
    return eig / norm, result ./ norm
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
            w = imag(iw[l])
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
    println("Getting Self Energy")
    Sigma = Self_Energy()
    println("Creating energy mesh")
    e = create_energy_mesh(band, iw, Sigma)

    println("Power Iteration to find Eigenvalue and symmetry")
    phi = rand(Float32, fnw, nx, ny, nz) 
    #phi = ones(fnw, nx, ny, nz) # Initialize 
    phi .= phi ./ sum(abs.(phi))^(0.5) # Normalize
    eig, phi = linearized_eliashberg_loop(phi, iw, V_rt, e, mesh)

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


end # module
