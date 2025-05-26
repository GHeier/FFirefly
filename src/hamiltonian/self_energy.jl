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
wc = cfg.bcs_cutoff_frequency
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


function zero_out_beyond_wc_e!(A, e)
    nw, nx, ny, nz = size(A)
    for j in 1:nx, k in 1:ny, l in 1:nz
        if abs(e[1, j, k, l]) > wc 
            A[:, j, k, l] .= 0
        end
    end
end

function get_G_V(iw, iv)
    fnw = length(iw)
    bnw = length(iv)
    G = Array{ComplexF32}(undef, fnw, nx, ny, nz)
    V = Array{ComplexF32}(undef, bnw, nx, ny, nz)
    e = Array{ComplexF32}(undef, 1, nx, ny, nz)
    println("Getting bands")
    band = Bands()
    println("Getting vertex")
    vertex = Vertex()
    println("Filling arrays")
    for iy in 1:ny, ix in 1:nx, iz in 1:nz
        kvec = get_kvec(ix - 1, iy - 1, iz - 1)
        e[1, ix, iy, iz] = band(1, kvec)
        for l in 1:fnw
            G[l, ix, iy, iz] = 1 / (iw[l] - band(1, kvec) )
        end
        for l in 1:bnw
            w = imag(iv[l])
            V[l, ix, iy, iz] = vertex(kvec, w)
        end
    end
    #zero_out_beyond_wc_e!(G, e)
    return G, V
end

function get_self_energy()
    is_bcs = false
    if interaction == "const"
        iw, iv = get_iw_iv_bcs()
        fnw, bnw = nw, nw
        is_bcs = true
        println("Taking const interaction")
    else
        println("Getting basis")
        basis = IR_Mesh()
        fnw, bnw = basis.fnw, basis.bnw
        iw, iv = get_iw_iv(basis)
    end

    G, V = get_G_V(iw, iv)
    println("Got arrays")
    if is_bcs
        p_fft = plan_fft(G)
        p_ifft = plan_ifft(G)

        G_rt = p_fft * G
        V_rt = p_fft * V
        V_rt .*= G_rt
        sigma = fftshift(p_ifft * V_rt) / (beta * nx * ny * nz)
        #sigma = ones(ComplexF32, nw, nx, ny, nz)
        N = Field_R(outdir * prefix * "_DOS.dat")
        println("Predicted BCS sigma = ", 1 * N(mu) * wc / 2)
    else
        G_rt = kw_to_rtau(G, 'F', basis)
        V_rt = kw_to_rtau(V, 'B', basis)
        temp = V_rt .* G_rt
        sigma = rtau_to_kw(temp, 'F', basis)
    end
    println("Max Sigma: ", maximum(real.(sigma)))
    println("Min Sigma: ", minimum(real.(sigma)))
    save(iw, sigma)
end

function save(iw_arr, sigma)
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
        fnw = length(iw_arr)
        for ix in 1:nx, iy in 1:ny, iz in 1:nz, l in 1:fnw
            kvec = get_kvec(ix - 1, iy - 1, iz - 1)
            iw = iw_arr[l]
            if dim == 3
                @printf(f, "%f %f %f %f %f %f\n", kvec[1], kvec[2], kvec[3], iw.im, sigma[l,ix,iy,iz].re, sigma[l,ix,iy,iz].im)
            elseif dim == 2
                @printf(f, "%f %f %f %f %f\n", kvec[1], kvec[2], iw.im, sigma[l,ix,iy,iz].re, sigma[l,ix,iy,iz].im)
            elseif dim == 1
                @printf(f, "%f %f %f %f\n", kvec[1], iw.im, sigma[l,ix,iy,iz].re, sigma[l,ix,iy,iz].im)
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


