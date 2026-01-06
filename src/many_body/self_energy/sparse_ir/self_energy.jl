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
using Interpolations

T = cfg.Temperature
beta = 1 / T
dim = cfg.dimension
q_mesh = cfg.q_mesh
nx, ny, nz = cfg.k_mesh
if dim == 2
    nz = 1
end
nw = cfg.w_pts
nbnd = cfg.nbnd
mu = cfg.fermi_energy
wc = cfg.cutoff_energy
BZ = cfg.brillouin_zone
prefix = cfg.prefix
outdir = cfg.outdir
interaction = cfg.interaction

function get_kvec(ix, iy, iz, nx, ny, nz)
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


function get_ek(band)
    e = Array{ComplexF32}(undef, 1, nx, ny, nz)
    for ix in 1:nx, iy in 1:ny, iz in 1:nz
        kvec = get_kvec(ix - 1, iy - 1, iz - 1, nx, ny, nz)
        e[1, ix, iy, iz] = band(1, kvec)
    end
    return e
end


function get_G_V(iw, iv, ek)
    fnw = length(iw)
    bnw = length(iv)
    G = Array{ComplexF32}(undef, fnw, nx, ny, nz)
    V = Array{ComplexF32}(undef, bnw, nx, ny, nz)
    println("Getting vertex")
    vertex = Vertex()
    println("Filling arrays")
    for iy in 1:ny, ix in 1:nx, iz in 1:nz
        kvec = get_kvec(ix - 1, iy - 1, iz - 1, nx, ny, nz)
        for l in 1:fnw
            G[l, ix, iy, iz] = 1 / (iw[l] - ek[1, ix, iy, iz])
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
    band = Bands()
    ek = get_ek(band)
    is_bcs = false
    if interaction == "const"
        iw, iv = get_iw_iv_bcs()
        fnw, bnw = nw, nw
        is_bcs = true
        println("Taking const interaction")
    else
        println("Getting basis")
        basis = IR_Mesh(ek)
        fnw, bnw = basis.fnw, basis.bnw
        iw, iv = get_iw_iv(basis)
    end

    G, V = get_G_V(iw, iv, ek)
    println("Got arrays")
    if is_bcs
        p_fft = plan_fft(G)
        p_ifft = plan_ifft(G)

        G_rt = p_fft * G
        V_rt = p_fft * V
        V_rt .*= G_rt
        sigma = fftshift(p_ifft * V_rt) / (beta * nx * ny * nz)
        #sigma = ones(ComplexF32, nw, nx, ny, nz)
        N = Field_R(outdir * prefix * "_DOS." + filetype) 
        println("Predicted BCS sigma = ", 1 * N(mu) * wc / 2)
    else
        G_rt = kw_to_rtau(G, 'F', basis)
        V_rt = kw_to_rtau(V, 'B', basis)
        temp = V_rt .* G_rt
        sigma = rtau_to_kw(temp, 'F', basis)
    end
    for i in 1:fnw
        sigma[i, :, :, :] = fftshift(sigma[i, :, :, :])
    end
    println("Max Sigma: ", maximum(real.(sigma)))
    println("Min Sigma: ", minimum(real.(sigma)))
    println("Interpolaing result")
    if nz > 1
        itp = interpolate((1:fnw, 1:nx, 1:ny, 1:nz), sigma, Gridded(Linear()))
    else
        sigma = repeat(sigma, 1, 1, 1, 2)
        itp = interpolate((1:fnw, 1:nx, 1:ny, 1:2), sigma, Gridded(Linear()))
    end
    println("Interpolated")
    save(iw, itp)
end

function save(iw_arr, sigma)
    println("Saving Self Energy")

    filepath = outdir * prefix * "_self_energy." + filetype
    max_SE = -Inf
    min_SE = Inf

    # Determine mesh size
    fnw = length(iw_arr)
    nx, ny, nz = q_mesh
    if dim == 2
        nz = 1
    end

    # Determine output header and formatter
    header, fmt = begin
        if dim == 3
            ("# x y z w Re(f) Im(f)\n", (k, iw, val) -> @sprintf("%f %f %f %f %f %f\n", k[1], k[2], k[3], iw.im, real(val), imag(val)))
        elseif dim == 2
            ("# x y w Re(f) Im(f)\n", (k, iw, val) -> @sprintf("%f %f %f %f %f\n", k[1], k[2], iw.im, real(val), imag(val)))
        elseif dim == 1
            ("# x w Re(f) Im(f)\n", (k, iw, val) -> @sprintf("%f %f %f %f\n", k[1], iw.im, real(val), imag(val)))
        else
            error("Invalid dimension: $dim")
        end
    end

    # Write data
    open(filepath, "w") do f
        print(f, header)
        for ix in 1:nx, iy in 1:ny, iz in 1:nz, l in 1:fnw
            kvec = get_kvec(ix - 1, iy - 1, iz - 1, nx, ny, nz)
            iw = iw_arr[l]
            val = sigma[l, ix, iy, iz]
            max_SE = max(real(val), max_SE)
            min_SE = min(real(val), min_SE)
            print(f, fmt(kvec, iw, val))
        end
    end

    println("Max Sigma Saved: ", max_SE)
    println("Min Sigma Saved: ", min_SE)
    println("Saved to ", filepath)
end

end


