module Eliashberg
println("Started Julia")
using Firefly

include("../objects/mesh.jl")
using .IRMesh

using SparseIR: Statistics, value, valueim, MatsubaraSampling64F, TauSampling64
using FastGaussQuadrature
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


function get_V(V, q, w1, w2)
    if abs(w1) > wc || abs(w2) > w2
        return 0
    end
    q = 
    return V(q, imag(w1 - w2))
end 


function fill_V_arr(iw, vertex, surfaces, band)
    num_surfaces = length(surfaces)
    fnw = length(iw)

    V = Matrix{Array{ComplexF32, 3}}(undef, num_surfaces, num_surfaces)

    for s1 in 1:num_surfaces
        i1 = length(surfaces[s1])
        for s2 in 1:num_surfaces
            i2 = length(surfaces[s2])
            V[s1, s2] = Array{ComplexF32}(undef, i1, i2, fnw)
            Vblock = V[s1, s2]

            @inbounds for i in 1:i1, j in 1:i2
                k1 = surfaces[s1][i]
                k2 = surfaces[s2][j]
                dA = k2.area
                vk = norm(Firefly.vk(k2.n, k2, band))

                for f in 1:fnw
                    Vblock[i, j, f] = get_V(vertex, k1 - k2, iw[f], 0) * dA / vk
                end
            end
        end
    end

    return V
end


function fill_F_G!(phi, Z, chi, iw, e, F, G)
    for l in 1:num_surfaces
        F[l] .= phi[l] ./ ((imag.(iw)' .* Z[l]).^2 .+ phi[l].^2 .+ (e' .- mu .+ chi[l]).^2)
        G[l] .= -(iw' .* Z[l] .+ e' .+ chi[l] .- mu) ./ (-(Z[l] .* iw').^2 .+ phi[l].^2 .+ (e' .+ chi[l] .- mu).^2)
    end
end


function compute_phi_sigma(phi, Z, chi, iw, e_4d, V)
    P1 = [similar(m) for m in phi]
    P2 = [similar(m) for m in phi]
    for i in 1:num_surfaces
        li = length(phi[i])
        e = e_4d[i]
        for j in 1:li
            phiw = phi[i][j]
            Zw = Z[i][j]
            chiw = chi[i][j]
            Fw = phiw ./ ((Zw .* iw).^2 .+ phiw .^ 2 .+ (e .- mu .+ chiw).^2)
            Gw = -(iw .* Zw .+ e .+ chiw .- mu) ./ (-(Zw .* iw).^2 .+ phiw.^2 .+ (e .+ chiw .- mu).^2)
            for a in 1:num_surfaces
                la = length(phi[a])
                for b in 1:la
                    Vw = V[i, a][j, b, :]
                    X1 = element_convolution_const(Vw, Fw)
                    X2 = element_convolution_const(Vw, Gw)
                    P1[i][j] += X1
                    P2[i][j] += X2
                end
            end
        end
    end
    return P1, P2
end


function fill_Z_chi!(sigma, iw, Z, chi)
    for l in 1:num_surfaces
        iw_2d = reshape(iw, 1, length(Z[l]))
        s = sigma[l]
        s_rev = s[:, end:-1:1]
        Z[l] = 1.0f0 .- 0.5f0 .* (s .- s_rev) ./ iw_2d
        chi[l] = 0.5 .* (s .+ s_rev)
    end
end
    

function compute!(V, phi, Z, chi, iw, sigma, e)
    phi, sigma = compute_phi_sigma(phi, Z, chi, iw, e, V)
    fill_Z_chi!(sigma, iw, Z, chi)
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


function initialize_ep_arr(surfaces, band, surface_vals, fnw)
    println("Initializing epsilon")
    e_4d = Vector{Matrix{Float32}}(undef, num_surfaces)
    for l in 1:num_surfaces
        i_len = length(surfaces[l])
        temp = ones(ComplexF32, i_len, fnw) * surface_vals[l]
        e_4d[l]  = repeat(temp, i_len, 1)  # broadcast-friendly
    end
    return e_4d
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


function initialize(iw, vertex, surfaces, band)
    println("Initializing Arrays")
    fnw = length(iw)
    phi = Vector{Matrix{ComplexF32}}(undef, num_surfaces)
    Z   = Vector{Matrix{ComplexF32}}(undef, num_surfaces)
    chi = Vector{Matrix{ComplexF32}}(undef, num_surfaces)
    sigma   = Vector{Matrix{ComplexF32}}(undef, num_surfaces)
    #F   = Vector{Matrix{ComplexF32}}(undef, num_surfaces)
    #G   = Vector{Matrix{ComplexF32}}(undef, num_surfaces)

    temp = 1e-3 ./ (imag.(iw).^2 .+ 1)  # fnw-length vector

    for l in 1:num_surfaces
        i_len = length(surfaces[l])
        phi[l]  = repeat(temp', i_len, 1)  # broadcast-friendly
        Z[l]    = ones(ComplexF32, i_len, fnw)
        chi[l]  = ones(ComplexF32, i_len, fnw)
        sigma[l]    = zeros(ComplexF32, i_len, fnw)
        #F[l]    = zeros(ComplexF32, i_len, fnw)
        #G[l]    = zeros(ComplexF32, i_len, fnw)

        #fill_F_G!(phi, Z, chi, iw, e_4d, F, G)
    end

    println("Initializing V")
    V = fill_V_arr(iw, vertex, surfaces, band)
    return phi, Z, chi, sigma, V
end


function converge(phi, Z, chi, sigma, V, iw, e_4d)
    println("Starting Convergence Loop")
    phierr = 0.0
    prev_phi = copy(phi)
    iterations = 1000
    for i in 1:iterations
        compute!(V, phi, Z, chi, iw, sigma, e_4d)

        phierr = minimum((maximum(abs.(real.(phi - prev_phi))), maximum(abs.(real.(phi + prev_phi)))))
        max_phi = round(maximum(abs.(phi)), digits=6)
        print("Iteration $i: MaxPhi = $max_phi              Error = $phierr           \r")

        if abs(phierr) < scf_tol || max_phi < 1e-10
            println("Converged after ", i, " iterations                                                                   ")
            break
        end
        prev_phi = copy(phi)
        #fill_F_G!(phi, Z, chi, iw, e_4d, F, G)
    end
    println("Error: ", phierr, "            ")
    return phi, Z, chi, sigma
end


function eliashberg_convsum()
    println("Beginning Eliashberg")
    scf_tol = 1e-4

    println("Getting Bands")
    band = Firefly.Bands()
    surfaces, surface_vals = get_surface(num_surfaces, band)

    println("Getting Vertex")
    vertex = Firefly.Vertex()
    println("Creating IRMesh")
    mesh = IR_Mesh()
    iw, iv = get_iw_iv(mesh)
    fnw, bnw = length(iw), length(iv)

    e_4d = initialize_ep_arr(surfaces, band, surface_vals, fnw)

    phi, Z, chi, sigma, V = initialize(iw, vertex, surfaces, band)

    phi, Z, chi, sigma = converge(phi, Z, chi, sigma, V, iw, e_4d)

    max_phi = maximum(abs.(phi))
    max_Z = maximum(abs.(Z))

    println("Max phi: ", max_phi)
    println("Max Z: ", max_Z)
    return
end


function eliashberg_node()
    eliashberg_convsum()
end

end # module

