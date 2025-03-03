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
prefix = cfg.prefix
outdir = cfg.outdir
dynamic = cfg.dynamic

function epsilon(kvec)
    if dim == 1
        return -2.0*cos(kvec[1])
    elseif dim == 2
        return -2.0*(cos(kvec[1])+cos(kvec[2]))
    elseif dim == 3
        return -2.0*(cos(kvec[1])+cos(kvec[2])+cos(kvec[3]))
    else
        println("Invalid dimension")
        exit(1)
    end
end

function get_kvec(ix, iy, iz)
    kvec = [ix / (nx - 1) - 0.5, iy / (ny - 1) - 0.5, iz / (nz - 1) - 0.5] 
    if dim < 3
        kvec[3] = 0.0
    elseif dim < 2
        kvec[2] = 0.0
    end
    kvec = BZ * kvec
    return kvec
end

function get_ckio_ir()
    println("Calculating IR response")
    basis = IR_Mesh(1e-15)
    ek = Array{ComplexF64,3}(undef, nx, ny, nz)
    for iy in 1:ny, ix in 1:nx, iz in 1:nz
        kvec = get_kvec(ix - 1, iy - 1, iz - 1)
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
    grio = k_to_r(gkio, basis)
    grit = wn_to_tau(basis, Fermionic(), grio)
    return grit
end

function ckio_calc(basis, grit)
    crit = Array{ComplexF64}(undef, basis.bntau, nx, ny, nz)
    for iy in 1:ny, ix in 1:nx, iz in 1:nz, it in 1:basis.bntau
        crit[it,ix,iy,iz] = grit[it,ix,iy,iz] * grit[basis.bntau-it+1,ix,iy,iz]
    end
    # Fourier transform
    ckit = r_to_k(crit, basis)
    ckio = tau_to_wn(basis, Bosonic(), ckit)
    return ckio
end

function save_ckio_ir(basis, ckio::Array{ComplexF64,4})
    println("Saving IR response")
    open(outdir * prefix * "_chi.dat", "w") do f
        if dim == 3 && dynamic == true
            @printf(f, "# x y z w Re(f) Im(f)\n")
        elseif dim == 2 && dynamic == true
            @printf(f, "# x y w Re(f) Im(f)\n")
        elseif dim == 1 && dynamic == true
            @printf(f, "# x w Re(f) Im(f)\n")
        elseif dim == 3 && dynamic == false
            @printf(f, "# x y z Re(f) Im(f)\n")
        elseif dim == 2 && dynamic == false
            @printf(f, "# x y Re(f) Im(f)\n")
        elseif dim == 1 && dynamic == false
            @printf(f, "# x Re(f) Im(f)\n")
        else
            println("Invalid dimension")
            exit(1)
        end
        for ix in 1:nx, iy in 1:ny, iz in 1:nz, iw in 1:basis.bnw
            kvec = get_kvec(ix - 1, iy - 1, iz - 1)
            iv = valueim(basis.IR_basis_set.smpl_wn_b.sampling_points[iw], beta)
            if dynamic == false && iv.im != 0.0
                continue
            end
            if dim == 3 && dynamic == true
                @printf(f, "%f %f %f %f %f %f\n", kvec[1], kvec[2], kvec[3], iv.im, ckio[iw,ix,iy,iz].re, ckio[iw,ix,iy,iz].im)
            elseif dim == 2 && dynamic == true
                @printf(f, "%f %f %f %f %f\n", kvec[1], kvec[2], iv.im, ckio[iw,ix,iy,iz].re, ckio[iw,ix,iy,iz].im)
            elseif dim == 1 && dynamic == true
                @printf(f, "%f %f %f %f\n", kvec[1], iv.im, ckio[iw,ix,iy,iz].re, ckio[iw,ix,iy,iz].im)
            elseif dim == 3 && dynamic == false
                @printf(f, "%f %f %f %f %f\n", kvec[1], kvec[2], kvec[3], ckio[iw,ix,iy,iz].re, ckio[iw,ix,iy,iz].im)
            elseif dim == 2 && dynamic == false
                @printf(f, "%f %f %f %f\n", kvec[1], kvec[2], ckio[iw,ix,iy,iz].re, ckio[iw,ix,iy,iz].im)
            elseif dim == 1 && dynamic == false
                @printf(f, "%f %f %f\n", kvec[1], ckio[iw,ix,iy,iz].re, ckio[iw,ix,iy,iz].im)
            else
                println("Invalid dimension")
                exit(1)
            end
        end
    end
end

end

