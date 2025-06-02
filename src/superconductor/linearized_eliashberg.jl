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
projs = cfg.projections

const beta = 1 / cfg.Temperature

const bcs_debug = true

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
    iters = 0
    while eig_err > 1e-5
        if !bcs_debug
            eig, phi = linearized_eliashberg(phi, iw, V_rt, e, mesh)
        else
            eig, phi = linearized_eliashberg_bcs(phi, iw, V_rt, e)
        end
        eig_err = abs(eig - prev_eig)
        prev_eig = eig
        print("Eig: ", round(eig, digits=6), " Error: ", round(eig_err, digits=6), "   \n")
        iters += 1
    end
    println("Iterations: ", iters)
    return eig, phi
end


function linearized_eliashberg(phi, iw, V_rt, e, mesh)
    F = -phi .* ( 1 ./ (iw .- e)) .* (1 ./ (-iw .- e))
    #F = -phi ./ ((Z .* imag.(iw)).^2 .+ e.^2)
    F_rt = kw_to_rtau(F, 'F', mesh)
    temp = V_rt .* F_rt
    newphi = rtau_to_kw(temp, 'F', mesh)

    norm = sum(newphi .* newphi)^(0.5)
    newphi .= newphi ./ norm

    result = phi .* newphi # Resultant eigenvector/value
    eig = sum(result) / sum(phi .* phi)
    return eig, newphi 
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

    eig = sum(result .* phi)
    norm = sum(phi .* phi)
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

    if bcs_debug
        N = Firefly.Field_R(outdir * prefix * "_DOS.dat")
        initial_expected = log(1.134 * wc * beta) * N(mu)
        println("Initial Expected eig: $initial_expected")
    end
    if projs != ""
        println("Finding eigs for $projs projections")
        projections = get_projections(fnw, nx, ny, nz)
        println("Projections Created")
        eigs, c_vals = linearized_eliashberg_projections(projections, iw, V_rt, e, mesh)
        println("Eigenvalue Projections Found")
        for i in eachindex(eigs)
            println("Eig for $(projs[i]) wave is $(eigs[i])")
        end
    end
    println("Power Iteration to find Eigenvalue and symmetry")
    phi = rand(Float32, fnw, nx, ny, nz) 
    eig, phi = linearized_eliashberg_loop(phi, iw, V_rt, e, mesh)
    @printf("Max Eig: %.6f\n", real(eig))

    #if !bcs_debug
    #    save_gap(phi, iw)
    #end
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


function get_projections(fnw, nx, ny, nz)
    num_projs = 2
    projections = Vector{Array{Float32}}(undef, num_projs)
    for i in 1:num_projs
        phi = Array{Float32}(undef, fnw, nx, ny, nz)
        for j in 1:nx, k in 1:ny, l in 1:nz
            kvec = get_kvec(nx, ny, nz)
            if projs[i] == 's'
                phi[:, j, k, l] .= 1.0
            elseif projs[i] == 'd'
                phi[:, j, k, l] .= cos(kvec[1]) - cos(kvec[2])
            end
        end
        phi ./= (sum(phi .* phi, dims=(2,3,4))) # Normalize in k-space
        projections[i] = phi
    end
    return projections
end
    


function linearized_eliashberg_projections(projections, iw, V_rt, e, mesh)
    num_projs = length(projections)
    eigs = zeros(num_projs)
    c_vals = Vector{Array{Float32}}(undef, num_projs)
    conv_thresh = 1e-4
    conv_iters = 1
    errs = ones(conv_iters)
    for i in 1:num_projs
        Delta_0 = projections[i]
        c_vals[i] = ones(size(Delta_0))
        prev_eig = Inf
        iter = 0
        err = true
        while err
            if bcs_debug
                zero_out_beyond_wc_e!(Delta_0, e)
            end

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
            w = Delta_1 ./ sum(Delta_1 .* projections[i], dims=(2,3,4))
            # Step 4
            eig = ( sum( (projections[i] .* w).^2 ) )^(0.5)
            # Step 5
            c_vals[i] .= w ./ eig
            eigs[i] = eig

            # Error checking
            errs[iter % conv_iters + 1] = abs(prev_eig - eig)
            err = any(e -> e > conv_thresh, errs)

            # Update for next iteration
            prev_eig = eig
            Delta_0 .= Delta_1
            iter += 1
            println("iter = $iter, eig = $eig")
        end
    end
    return eigs, c_vals
end

end # module
