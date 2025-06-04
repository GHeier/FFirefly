module response_ir

using Firefly
cfg = Firefly.Config
#Firefly.load_config!("/home/g/Research/ffirefly/build/bin/input.cfg")

include("../objects/mesh.jl")
using .IRMesh
using SparseIR
import SparseIR: Statistics, value, valueim
using PyCall
using Printf
using Base.Threads
using Interpolations

T = cfg.Temperature
beta = 1 / T
dim = cfg.dimension
nx, ny, nz = cfg.k_mesh
if nx < cfg.q_mesh[1] || ny < cfg.q_mesh[2] || nz < cfg.q_mesh[3]
    println("Q Mesh cannot be larger than K Mesh.\nQ Mesh is the mesh where the data is saved")
    exit()
end
if dim == 2
    nz = 1
end
n = cfg.w_pts
mu = cfg.fermi_energy
BZ = cfg.brillouin_zone
prefix = cfg.prefix
outdir = cfg.outdir
dynamic = cfg.dynamic

function get_kvec(ix, iy, iz, nx, ny, nz)
    kvec = [ix / (nx) - 0.0, iy / (ny) - 0.0, iz / (nz) - 0.0] 
    if dim < 3
        kvec[3] = 0.0
    elseif dim < 2
        kvec[2] = 0.0
    end
    kvec = BZ * kvec
    return kvec
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


function get_ckio_ir()
    println("Calculating IR response")
    println("Filling bands array")
    band = Bands()
    ek = Array{Float64}(undef, 1, nx, ny, nz)
    for ix in 1:nx, iy in 1:ny, iz in 1:nz
        kvec = get_kvec(ix - 1, iy - 1, iz - 1, nx, ny, nz)
        ek[1, ix, iy, iz] = band(1, kvec)
    end
    println("Getting basis")
    basis = IR_Mesh(ek, 1e-15)
    iw, iv = get_iw_iv(basis)
    iw = reshape(iw, :, 1, 1, 1)


    println("Calculating G(k,iw)")
    gkio = gkio_calc(basis, ek, mu, iw)
    println("Calculating G(r,tau)")
    grit = grit_calc(basis, gkio)
    println("Calculating X(k,iv)")
    ckio = ckio_calc(basis, grit)
    println(ckio[Int(basis.bnw / 2 + 0.5), 1, 1, 1])
    println("Max χ = ", maximum(real.(ckio)))
    println("Min χ: ", minimum(real.(ckio)))
    println("Interpolaing result")
    if nz > 1
        itp = interpolate((1:basis.bnw, 1:nx, 1:ny, 1:nz), ckio, Gridded(Linear()))
    else
        ckio = repeat(ckio, 1, 1, 1, 2)
        itp = interpolate((1:basis.bnw, 1:nx, 1:ny, 1:2), ckio, Gridded(Linear()))
    end
    println("Interpolated")
    save_ckio_ir(basis, ckio)
end

function gkio_calc(basis, ek::Array{Float64,4}, mu::Float64, iw)
    gkio = 1.0 ./ (iw .+ mu .- ek)
    return gkio
end

function grit_calc(basis, gkio)
    grio = k_to_r(gkio, basis)
    grit = wn_to_tau(basis, Fermionic(), grio)
    return grit
end

function ckio_calc(basis, grit)
    crit = Array{ComplexF64}(undef, basis.bntau, nx, ny, nz)
    #crit .= grit .* grit
    crit .= grit .* reverse(grit, dims=1)
    # Fourier transform
    ckit = r_to_k(crit, basis)
    ckio = tau_to_wn(basis, Bosonic(), ckit)
    return ckio
end

function save_ckio_ir(basis, ckio)
    println("Saving IR response")
    max_X = -1000
    min_X = 1000
    nx, ny, nz = cfg.q_mesh
    if (dim == 2)
        nz = 1
    end
    open(outdir * prefix * "_chi.dat", "w") do f
        header, fmt = begin
            if dim == 3
                dynamic ? ("# x y z w Re(f) Im(f)\n", (k, iv, v) -> @sprintf("%f %f %f %f %f %f\n", k[1], k[2], k[3], iv.im, real(v), imag(v))) : ("# x y z Re(f) Im(f)\n",        (k, iv, v) -> @sprintf("%f %f %f %f %f\n",   k[1], k[2], k[3], real(v), imag(v)))
            elseif dim == 2
                dynamic ? ("# x y w Re(f) Im(f)\n", (k, iv, v) -> @sprintf("%f %f %f %f %f\n",   k[1], k[2], iv.im, real(v), imag(v))) : ("# x y Re(f) Im(f)\n",    (k, iv, v) -> @sprintf("%f %f %f %f\n",     k[1], k[2], real(v), imag(v)))
            elseif dim == 1
                dynamic ? ("# x w Re(f) Im(f)\n",   (k, iv, v) -> @sprintf("%f %f %f %f\n",     k[1], iv.im, real(v), imag(v))) : ("# x Re(f) Im(f)\n",      (k, iv, v) -> @sprintf("%f %f %f\n",       k[1], real(v), imag(v)))
            else
                error("Invalid dimension: $dim")
            end
        end

       print(f, header)

       for ix in 1:nx, iy in 1:ny, iz in 1:nz, iw in 1:basis.bnw
            kvec = get_kvec(ix - 1, iy - 1, iz - 1, nx, ny, nz)
            iv = valueim(basis.IR_basis_set.smpl_wn_b.sampling_points[iw], beta)
            if !dynamic && iv.im != 0.0
                continue
            end
            val = ckio[iw, ix, iy, iz]
            max_X = max(real(val), max_X)
            min_X = min(real(val), min_X)
            print(f, fmt(kvec, iv, val))
        end
    end
    println("Max χ Saved: ", max_X)
    println("Min χ Saved: ", min_X)
    println("Saved to ", outdir * prefix * "_chi.dat")
end

end

