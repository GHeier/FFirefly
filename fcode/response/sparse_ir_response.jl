module response_ir

include("../objects/mesh.jl")
using .IRMesh
using SparseIR
import SparseIR: Statistics, value, valueim
using PyCall
using Printf
fcode = pyimport("fcode")
cfg = fcode.config

T = cfg.Temperature
beta = 1 / T
dim = cfg.dimension
nx, ny, nz = cfg.k_mesh
if dim == 2
    nz = 1
end
n = cfg.w_pts
mu = cfg.fermi_energy
BZ = cfg.brillouin_zone

function epsilon(kvec)
    return -2.0*(cos(kvec[1])+cos(kvec[2]))
end

function get_ckio_ir()
    println("Calculating IR response")
    basis = IR_Mesh()
    ek = Array{ComplexF64,3}(undef, nx, ny, nz)
    for iy in 1:ny, ix in 1:nx, iz in 1:nz
        kvec = [ix / nx, iy / ny, iz / nz]
        kvec = BZ * kvec
        ek[ix, iy, iz] = epsilon(kvec)
    end

    gkio = gkio_calc(basis, ek, mu)
    grit = grit_calc(basis, gkio)
    ckio = ckio_calc(basis, grit)
    save_ckio_ir(basis, ckio)
end

function gkio_calc(basis, ek::Array{ComplexF64,3}, mu::Float64)
    gkio = Array{ComplexF64}(undef, basis.fnw, nx, ny, nz)
    for iy in 1:ny, ix in 1:nx, iw in 1:basis.fnw, iz in 1:nz
        iv::ComplexF64 = valueim(basis.IR_basis_set.smpl_wn_f.sampling_points[iw], beta)
        gkio[iw,ix,iy,iz] = 1.0/(iv - ek[ix,iy,iz] + mu)
    end
    return gkio
end

function grit_calc(basis, gkio)
    grio = k_to_r(basis, gkio)
    grit = wn_to_tau(basis, Fermionic(), grio)
    return grit
end

function ckio_calc(basis, grit)
    crit = Array{ComplexF64}(undef, basis.bntau, nx, ny, nz)
    for iy in 1:ny, ix in 1:nx, iz in 1:nz, it in 1:basis.bntau
        crit[it,ix,iy,iz] = grit[it,ix,iy,iz] * grit[basis.bntau-it+1,ix,iy,iz]
    end
    # Fourier transform
    ckit = r_to_k(basis, crit)
    ckio = tau_to_wn(basis, Bosonic(), ckit)
    return ckio
end

function save_ckio_ir(basis, ckio::Array{ComplexF64,4})
    println("Saving IR response")
    open("ckio_ir.dat", "w") do f
        for iy in 1:ny, ix in 1:nx, iz in 1:nz, iw in 1:basis.fnw
            kvec = [ix / nx, iy / ny, iz / nz]
            kvec = BZ * kvec
            iv = valueim(basis.IR_basis_set.smpl_wn_b.sampling_points[iw], beta)
            if dim == 3
                @printf(f, "%f %f %f %f %f %f\n", kvec[1], kvec[2], kvec[3], iv.im, ckio[iw,ix,iy,iz].re, ckio[iw,ix,iy,iz].im)
            elseif dim == 2
                @printf(f, "%f %f %f %f %f\n", kvec[1], kvec[2], iv.im, ckio[iw,ix,iy,iz].re, ckio[iw,ix,iy,iz].im)
            elseif dim == 1
                @printf(f, "%f %f %f %f\n", kvec[1], iv.im, ckio[iw,ix,iy,iz].re, ckio[iw,ix,iy,iz].im)
            else
                println("Invalid dimension")
                exit(1)
            end
        end
    end
end

end

using .response_ir
response_ir.get_ckio_ir()
