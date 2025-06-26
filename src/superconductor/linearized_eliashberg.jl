module Linearized_Eliashberg
println("Started Julia")
using Firefly

using LinearAlgebra: norm, dot
include("../objects/mesh.jl")
using .IRMesh

using SparseIR
import SparseIR: Statistics, value, valueim, MatsubaraSampling64F, TauSampling64
using Roots, Printf, DelimitedFiles
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
projs = cfg.projections

const verbosity = cfg.verbosity
const filetype = cfg.filetype

const beta = 1 / cfg.Temperature

const bcs_debug = false
if bcs_debug
    printstyled("BCS Debug Session enabled\n"; color = :blue)
    if nw % 2 != 0
        println("Remember to make nw even")
        exit()
    end
end

function get_kvec(i, j, k)
    temp = BZ * [i / nx - 0.0, j / ny - 0.0, k / nz - 0.0]
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
    for i in eachindex(iw)
        if abs(imag(iw[i])) > wc 
            A[i, :, :, :] .= 0
        end
    end
end

function fill_V_rt(vertex, iv, mesh)
    bnw = length(iv)
    V = Array{ComplexF32}(undef, bnw, nx, ny, nz)
    for i in 1:nx, j in 1:ny, k in 1:nz, l in 1:bnw
        q = get_kvec(i - 1, j - 1, k - 1)
        w = imag(iv[l])
        V[l, i, j, k] = 0.5 * (vertex(q, w) + vertex(-q, w))
    end
    if bcs_debug
        return fft(V)
    end
    return kw_to_rtau(V, 'B', mesh)
end

function create_energy_mesh(band, iw, Sigma, with_sigma=true)
    nw = length(iw)
    e_arr = zeros(ComplexF32, nw, nx, ny, nz)
    for i in 1:nx, j in 1:ny, k in 1:nz, l in 1:nw
        kvec = get_kvec(i - 1, j - 1, k - 1)
        if with_sigma
            e_arr[l, i, j, k] = band(1, kvec) - mu + Sigma(kvec, imag(iw[l]))
        else
            e_arr[l, i, j, k] = band(1, kvec) - mu
        end
    end
    return e_arr
end

function deflate_v!(v, eigvecs)
    for phi_i in eigvecs
        proj = (dot(vec(phi_i), vec(v)) / dot(vec(phi_i), vec(phi_i))) * phi_i
        v .-= proj
    end
    nv = norm(vec(v))
    if nv > 0
        v ./= nv
    else
        error("Deflation resulted in zero vector")
    end
end

function power_iteration(iw, V_rt, e, mesh = 0)
    println("Beginning Power Iteration")
    deflate = Vector{Array{ComplexF32, 4}}()
    eig, iters, max_iters = -1, 0, 5
    while real(eig) < 0 && iters < max_iters
        eig, phi = linearized_eliashberg_loop(iw, V_rt, e, deflate, mesh)
        println("eig$iters = $eig")
        save!(outdir * prefix * "_gap_$iters.h5", phi, iw)
        if eig != 0
            push!(deflate, conj.(phi))
        end
        iters += 1
    end
    return eig
end

function save!(filename, phi, iw)
    mesh = [nx, ny, nz]
    BZ = cfg.brillouin_zone
    if dim == 2
        mesh = mesh[1:end-1]
        BZ = BZ[1:end-1, 1:end-1]
    end
    save_field!(filename, phi, BZ, mesh, imag.(iw))
end


function linearized_eliashberg_loop(iw, V_rt, e, deflates, mesh = 0)
    phinw = length(iw)
    phi = rand(ComplexF32, phinw, nx, ny, nz)
    for iy in 1:ny, ix in 1:nx, iw in 1:phinw
        kx, ky = get_kvec(ix, iy, 1)
        #kx::Float64 = (2*π*(ix-1))/nx
        #ky::Float64 = (2*π*(iy-1))/ny
        #phi[iw,ix,iy,1] = cos(kx) - cos(ky)
    end
    eig, prev_eig, eig_err = 0, 0, 1
    max_iters = 50
    iters = 0
    while eig_err > 1e-4 && iters < max_iters
        if !bcs_debug
            eig, phi = linearized_eliashberg(phi, iw, V_rt, e, mesh)
        else
            eig, phi = linearized_eliashberg_bcs(phi, iw, V_rt, e)
        end
        deflate_v!(phi, deflates)
        eig_err = abs(eig - prev_eig)
        prev_eig = eig
        if verbosity == "high"
            print("Eig: ", round(real(eig), digits=6), " Error: ", round(eig_err, digits=6), "       \n")
        else
            print("Eig: ", round(real(eig), digits=6), " Error: ", round(eig_err, digits=6), "       \r")
        end
        iters += 1
    end
    println("Iterations: $iters                                 ")
    return eig, phi
end


function convolution(A, B, mesh)
    if bcs_debug
        return fftshift(ifft(fft(A) .* B)) / (beta * nk)
    else
        A_rt = kw_to_rtau(A, 'F', mesh)
        A_rt .= B .* A_rt
        temp = rtau_to_kw(A_rt, 'F', mesh)
        return temp
    end
end


function linearized_eliashberg(phi, iw, V_rt, e, mesh)
    F = -phi .* ( 1 ./ (iw .- e)) .* (1 ./ (-iw .- e))
    result = convolution(F, V_rt, mesh)

    eig = sum(real.(conj.(result) .* phi))
    norm = sum(real.(phi .* conj.(phi)))
    return eig / norm, result ./ norm^(0.5)
end


function linearized_eliashberg_bcs(phi, iw, V_rt, e)
    #zero_out_beyond_wc_e!(phi, e)
    #zero_out_beyond_wc_iw!(phi, iw)
    F = -phi .* ( 1 ./ (iw .- e)) .* (1 ./ (-iw .- e))
    result = convolution(F, V_rt, 0)
    #zero_out_beyond_wc_e!(phi, e)
    #zero_out_beyond_wc_iw!(result, iw)

    eig = sum(real.(conj.(result) .* phi))
    norm = sum(real.(phi .* conj.(phi)))
    return eig / norm, result ./ norm^(0.5)
end


function save_gap(phi, iw)
    println("Saving Gap")
    nw = length(iw)
    filename = outdir * prefix * "_gap." * filetype
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

    println("Getting Bands")
    band = Firefly.Bands()
    println("Creating energy mesh")
    e = create_energy_mesh(band, 0, 0, false)

    if !bcs_debug
        println("Creating Mesh")
        mesh = IR_Mesh(real.(e))
        iw, iv = get_iw_iv(mesh)
        fnw, bnw = length(iw), length(iv)
    else 
        mesh = 0
        iw, iv = get_iw_iv_bcs()
        fnw, bnw = nw, nw
        println("BCS DEBUG: Found iw/iv")
    end

    println("Getting Self Energy")
    Sigma = Self_Energy()
    e = create_energy_mesh(band, iw, Sigma, true)

    if !bcs_debug
        println("Getting Vertex")
        vertex = Firefly.Vertex()
        println("Fourier Transforming Vertex")
        V_rt = fill_V_rt(vertex, iv, mesh)
        save!("V_rt_si.h5", V_rt, iv)
    else
        #V = -1.0 .* ones(Complex{Float32}, bnw, nx, ny, nz)
        #V_rt = fft(V)
        vertex = Firefly.Vertex()
        V_rt = fill_V_rt(vertex, iv, 0)
        println("BCS DEBUG: FFT'd Vertex")
        save!("V_rt_mf.h5", V_rt, iv)
    end

    if bcs_debug
        N = Firefly.Field_R(outdir * prefix * "_DOS." * filetype)
        initial_expected = log(1.134 * wc * beta) * N(mu)
        println("Initial Expected eig: $initial_expected")
    end
    if projs != ""
        println("Finding eigs for $projs projections")
        projections = get_projections(fnw, nx, ny, nz, length(projs))
        println("Projections Created")
        eigs, c_vals = linearized_eliashberg_projections(projections, iw, V_rt, e, mesh)
        println("Eigenvalue Projections Found")
        for i in eachindex(eigs)
            @printf("Eig for %c wave is %.4f\n", projs[i], eigs[i])
        end
    end
    println("Power Iteration to find Eigenvalue and symmetry")
    eig = power_iteration(iw, V_rt, e, mesh)
    #eig, phi = linearized_eliashberg_loop(iw, V_rt, e, eigs, mesh)
    @printf("Max Eig: %.4f\n", real(eig))
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


function get_projections(fnw, nx, ny, nz, num_projs)
    projections = Vector{Array{Float32}}(undef, num_projs)
    for i in 1:num_projs
        phi = Array{Float32}(undef, fnw, nx, ny, nz)
        for j in 1:nx, k in 1:ny, l in 1:nz
            kvec = get_kvec(j, k, l)
            if projs[i] == 's'
                phi[:, j, k, l] .= 1.0
            elseif projs[i] == 'd'
                phi[:, j, k, l] .= cos(kvec[1]) - cos(kvec[2])
            end
        end
        #phi ./= (sum(phi .* phi)^(0.5))
        projections[i] = phi
    end
    return projections
end



function linearized_eliashberg_projections(projections, iw, V_rt, e, mesh)
    nw = length(iw)
    num_projs = length(projections)
    eigs = zeros(num_projs)
    c_vals = Vector{Array{Float32}}(undef, num_projs)
    conv_thresh = 1e-4
    conv_iters = 2
    errs = ones(conv_iters)
    for i in 1:num_projs
        Delta_0 = projections[i]
        c_vals[i] = ones(size(Delta_0))
        iter = 0
        err = true
        while err
            if bcs_debug
                zero_out_beyond_wc_e!(Delta_0, e)
            end
            Delta_0 ./= (sum(Delta_0 .* Delta_0)^(0.5)) # Normalize Before calculation

            # Step 1
            F = -( 1 ./ (iw .- e)) .* (1 ./ (-iw .- e)) .* Delta_0 .* c_vals[i]
            # Step 2
            if !bcs_debug
                F_rt = kw_to_rtau(F, 'F', mesh)
                F_rt .= V_rt .* F_rt
                Delta_1 = real.(rtau_to_kw(F_rt, 'F', mesh))
            else
                F_rt = fft(F)
                F_rt .= V_rt .* F_rt
                Delta_1 = real.(fftshift(ifft(F_rt)) ./ (nk * beta))
                zero_out_beyond_wc_e!(Delta_1, e)
            end

            # Step 3
            eig = sum(Delta_1 .* Delta_0)
            # Step 4
            c_n = sum(Delta_1 .* Delta_0, dims=(2, 3, 4)) ./ eig .* nw
            printv("max c_n: $(maximum(c_n))")

            # Error checking
            errs[iter % conv_iters + 1] = abs(eigs[i] - eig)
            err = any(e -> e > conv_thresh, errs)

            # Update for next iteration
            c_vals[i] = c_n
            eigs[i] = eig
            Delta_0 .= Delta_1
            iter += 1
            printv("iter = $iter, eig = $eig")
        end
    end
    return eigs, c_vals
end

end # module
