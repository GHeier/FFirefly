module Linearized_Eliashberg
println("Started Julia")
using Firefly

include("../objects/mesh.jl")
using .IRMesh

using SparseIR: Statistics, value, valueim, MatsubaraSampling64F, TauSampling64
using Roots
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

function get_kvec(i, j, k)
    temp = BZ * [i / nx - 0.5, j / ny - 0.5, k / nz - 0.5]
    return temp[1:dim]
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
        F = iw .* Z ./ ( (iw .* Z) .^2 .+ e.^2)
        F_rt = kw_to_rtau(F, 'F', mesh)

        temp = V_rt .* F_rt
        result = rtau_to_kw(temp, 'F', mesh)

        Z .= 1 .+ result ./ iw

        Z_max_f = maximum(abs.(real.(Z)))
        Z_err = abs(Z_max_f - Z_max_i)
        Z_max_i = Z_max_f
        println("Z_max: ", Z_max_f)
    end
end


function linearized_eliashberg(phi, Z, iw, V_rt, e, mesh)
    F = -phi ./ ((Z .* iw).^2 .+ e.^2)
    F_rt = kw_to_rtau(F, 'F', mesh)
    temp = V_rt .* F_rt
    result = rtau_to_kw(temp, 'F', mesh)
    return maximum(abs.(real.(result))) / maximum(abs.(real.(phi))), result
end


function save_gap(phi, iw)
    nw = length(iw)
    filename = outdir * prefix * "_gap.dat"
    open(filename, "w") do io
        if dim == 3
            println(io, "kx\tky\tkz\tw\tRe(f)\tIm(f)")
        else 
            println(io, "kx\tky\tw\tRe(f)\tIm(f)")
        end


        for l in 1:nw, i in 1:nx, j in 1:ny, k in 1:nz
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
    println("Creating Mesh")
    mesh = IR_Mesh()
    iw, iv = get_iw_iv(mesh)
    fnw, bnw = length(iw), length(iv)

    println("Getting Vertex")
    vertex = Firefly.Vertex()
    println("Fourier Transforming Vertex")
    V_rt = fill_V_rt(vertex, iv, mesh)

    println("Getting Bands")
    band = Firefly.Bands()
    println("Creating energy mesh")
    e_4d = create_energy_mesh(band)

    println("Calculating Z")
    Z = ones(Complex{Float32}, fnw, nx, ny, nz)
    Z_convergence!(Z, V_rt, iw, e_4d, mesh)

    println("Calculating Eigenvalue and symmetry")
    phi = ones(Complex{Float32}, fnw, nx, ny, nz)
    eig, phi = linearized_eliashberg(phi, Z, iw, V_rt, e_4d, mesh)
    println("Eig: ", eig)
    save_gap(phi, iw)
end


end # module
