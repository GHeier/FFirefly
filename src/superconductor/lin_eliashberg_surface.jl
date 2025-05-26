module Linearized_Eliashberg_Surface
println("Started Julia")
using Firefly

include("../objects/mesh.jl")
using .IRMesh

using SparseIR
import SparseIR: Statistics, value, valueim, MatsubaraSampling64F, TauSampling64
using Roots, Printf, DelimitedFiles
using FastGaussQuadrature
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

const num_surfaces = 7


function get_kvec(i, j, k)
    temp = BZ * [i / nx - 0.5, j / ny - 0.5, k / nz - 0.5]
    return temp[1:dim]
end


function get_V(V, q, w)
    l = length(w)
    return [V(q, imag(w[i])) for i in 1:l]
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


function initialize_ep_phi_arr(surfaces, surface_vals, fnw)
    println("Initializing epsilon")
    e_4d = Vector{Float32}(undef, num_surfaces)
    phi = Vector{Matrix{ComplexF32}}(undef, num_surfaces)
    for l in 1:num_surfaces
        i_len = length(surfaces[l])
        e_4d[l] = surface_vals[l]
        phi[l] = ones(ComplexF32, i_len, fnw)
    end
    return e_4d, phi
end


function fill_F_G!(phi, Z, chi, iw, e, F, G)
    for l in 1:num_surfaces
        F[l] .= phi[l] ./ ((imag.(iw)' .* Z[l]).^2 .+ phi[l].^2 .+ (e' .- mu .+ chi[l]).^2)
        G[l] .= -(iw' .* Z[l] .+ e' .+ chi[l] .- mu) ./ (-(Z[l] .* iw').^2 .+ phi[l].^2 .+ (e' .+ chi[l] .- mu).^2)
    end
end


function element_convolution_F(A, B, mesh)
    A_t = wn_to_tau(mesh, Bosonic(), A)
    B_t = wn_to_tau(mesh, Fermionic(), B)
    temp = A_t .* B_t
    return tau_to_wn(mesh, Fermionic(), temp)
end


function fill_V_arr(iv, vertex, surfaces, band)
    num_surfaces = length(surfaces)
    bnw = length(iv)

    V_arr = Matrix{Array{ComplexF32, 3}}(undef, num_surfaces, num_surfaces)

    for s1 in 1:num_surfaces
        i1 = length(surfaces[s1])
        for s2 in 1:num_surfaces
            i2 = length(surfaces[s2])
            Vblock = Array{ComplexF32}(undef, i1, i2, bnw)

            @inbounds for i in 1:i1, j in 1:i2
                k1 = surfaces[s1][i]
                k2 = surfaces[s2][j]
                dA = k2.area
                vk_mag = Firefly.norm(vk(k2.n, k2, band))

                Vblock[i, j, :] = get_V(vertex, k1 - k2, iv) * dA / vk_mag
            end
            V_arr[s1, s2] = Vblock
        end
    end

    return V_arr
end


function compute_phi(phi, iw, e_4d, V, mesh)
    P1 = [similar(m) for m in phi]
    num_k = 0
    for i in 1:num_surfaces
        li = size(phi[i], 1)
        num_k += li
        e = e_4d[i]
        for j in 1:li
            phiw = phi[i][j, :]
            Fw = -phiw .* ( 1 ./ (iw .- e)) .* (1 ./ (-iw .- e))
            b = false
            for x in Fw
                b = isnan(x)
                if b 
                    break 
                end
            end
            if b
                println("Fw NAN: ", Fw)
            end
            for a in 1:num_surfaces
                la = size(phi[a], 1)
                for b in 1:la
                    Vw = V[i, a][j, b, :]
                    X1 = element_convolution_F(Vw, Fw, mesh)
                    b = false
                    for x in Vw
                        b = isnan(x)
                        if b 
                            break 
                        end
                    end
                    if b
                        #println("Vw NAN: ", Vw)
                    end
                    P1[i][j, :] += X1
                end
            end
        end
    end
    surf_eigs = [sum(P1[i] .* phi[i]) for i in 1:num_surfaces]
    eig = sum(surf_eigs) / (num_k * beta)
    return eig, P1
end


function get_surface(num_surfaces, band)
    println("Calculating Surfaces")
    surface_vals, surface_weights = gausslegendre(num_surfaces)
    surfaces = Vector{Vector{Vec}}(undef, num_surfaces)
    for i in 1:num_surfaces, n in 1:nbnd
        function band_func(k)
            return band(n, k)
        end
        surfaces[i] = Firefly.Surface(band_func, wc * surface_vals[i]).faces
    end
    return surfaces, surface_vals * wc
end


function linearized_eliashberg_loop(phi, iw, V_rt, e, mesh)
    eig, prev_eig, eig_err = 0, 0, 1
    while eig_err > 1e-5
        println("1")
        eig, phi = compute_phi(phi, iw, e, V_rt, mesh)
        println("1")
        println(eig)
        eig_err = abs(eig - prev_eig)
        prev_eig = eig
        println("Eig: ", round(eig, digits=6), " Error: ", round(eig_err, digits=6), "   \r")
        return eig, phi
    end
    println()
    return eig, phi
end


function eigenvalue_computation()
    println("Started Julia")
    println("Beginning Computing eigenvalues")
    println("Creating Mesh")
    mesh = IR_Mesh()
    iw, iv = get_iw_iv(mesh)
    fnw, bnw = length(iw), length(iv)

    println("Getting Bands")
    band = Firefly.Bands()
    surfaces, surface_vals = get_surface(num_surfaces, band)
    println("Getting Vertex")
    vertex = Firefly.Vertex()
    println("Creating V(k1, k2, tau)")
    V_rt = fill_V_arr(iv, vertex, surfaces, band)


    println("Creating energy mesh")
    e_4d, phi = initialize_ep_phi_arr(surfaces, surface_vals, fnw)

    println("Iteration to find Eigenvalue and symmetry")
    eig, phi = linearized_eliashberg_loop(phi, iw, V_rt, e_4d, mesh)
    @printf("Final Eig: %.6f\n", real(eig))
end

end # module
