module Self_Energy

using Firefly
cfg = Firefly.Config
#Firefly.load_config!("/home/g/Research/ffirefly/build/bin/input.cfg")

include("../objects/mesh.jl")
using .IRMesh
using SparseIR
import SparseIR: Statistics, value, valueim
using FFTW
using PyCall
using Printf
using Base.Threads

T = cfg.Temperature
beta = 1 / T
dim = cfg.dimension
nx, ny, nz = cfg.k_mesh
if dim == 2
    nz = 1
end
nw = cfg.w_pts
nbnd = cfg.nbnd
mu = cfg.fermi_energy
BZ = cfg.brillouin_zone
prefix = cfg.prefix
outdir = cfg.outdir
interaction = cfg.interaction

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

function get_self_energy()
    println("Getting vertex")
    vertex = Vertex()
    k1 = [0.1, 0.2, 0.3]
    k2 = [0.2, 0.3, 0.4]
    is_bcs = false
    if abs(vertex(k1) - vertex(k2)) < 1e-4 && interaction == "const"
        iw, iv = get_iw_iv_bcs()
        fnw, bnw = nw, nw
        is_bcs = true
    else
        println("Getting basis")
        basis = IR_Mesh(1e-15)
        fnw, bnw = basis.fnw, basis.bnw
        iw, iv = get_iw_iv(basis)
    end

    println("Creating G and V arrays")
    G = Array{ComplexF64}(undef, fnw, nx, ny, nz)
    V = Array{ComplexF64}(undef, bnw, nx, ny, nz)
    println("Getting bands")
    band = Bands()
    println("Fill arrays")
    for iy in 1:ny, ix in 1:nx, iz in 1:nz
        kvec = get_kvec(ix - 1, iy - 1, iz - 1)
        for l in 1:fnw
            w = imag(iw[l])
            G[l, ix, iy, iz] = band(1, kvec) 
        end
        for l in 1:bnw
            w = imag(iv[l])
            V[l, ix, iy, iz] = vertex(kvec, w)
        end
    end
    if is_bcs
        G_rt = fft(G)
        V_rt = fft(V)
        temp = V_rt .* G_rt
        sigma = fftshift(ifft(temp)) / (beta * nx * ny * nz)
    else
        G_rt = kw_to_rtau(G, "F", basis)
        V_rt = kw_to_rtau(V, "B", basis)
        temp = V_rt .* G_rt
        sigma = rtau_to_kw(temp, "B", basis)
    end
    println("Max Sigma: ", maximum(real.(sigma)))
    println("Min Sigma: ", minimum(real.(sigma)))
    #save(iw, sigma)
end

function save(iv_arr, sigma)
    println("Saving Self Energy")
    open(outdir * prefix * "_self_energy.dat", "w") do f
        if dim == 3
            @printf(f, "# x y z w Re(f) Im(f)\n")
        elseif dim == 2
            @printf(f, "# x y w Re(f) Im(f)\n")
        elseif dim == 1
            @printf(f, "# x w Re(f) Im(f)\n")
        else
            println("Invalid dimension")
            exit(1)
        end
        bnw = length(iv_arr)
        for ix in 1:nx, iy in 1:ny, iz in 1:nz, l in 1:bnw
            kvec = get_kvec(ix - 1, iy - 1, iz - 1)
            iv = iv_arr[l]
            if dim == 3
                @printf(f, "%f %f %f %f %f %f\n", kvec[1], kvec[2], kvec[3], iv.im, sigma[l,ix,iy,iz].re, sigma[l,ix,iy,iz].im)
            elseif dim == 2
                @printf(f, "%f %f %f %f %f\n", kvec[1], kvec[2], iv.im, sigma[l,ix,iy,iz].re, sigma[l,ix,iy,iz].im)
            elseif dim == 1
                @printf(f, "%f %f %f %f\n", kvec[1], iv.im, sigma[l,ix,iy,iz].re, sigma[l,ix,iy,iz].im)
            elseif dim == 3
                @printf(f, "%f %f %f %f %f\n", kvec[1], kvec[2], kvec[3], sigma[l,ix,iy,iz].re, sigma[l,ix,iy,iz].im)
            elseif dim == 2
                @printf(f, "%f %f %f %f\n", kvec[1], kvec[2], sigma[l,ix,iy,iz].re, sigma[l,ix,iy,iz].im)
            elseif dim == 1
                @printf(f, "%f %f %f\n", kvec[1], sigma[l,ix,iy,iz].re, sigma[l,ix,iy,iz].im)
            else
                println("Invalid dimension")
                exit(1)
            end
        end
    end
    println("Saved to ", outdir * prefix * "_self_energy.dat")
end

end


