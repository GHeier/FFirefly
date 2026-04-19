using SparseIR
import SparseIR: Statistics, value, valueim
using FFTW
using LinearAlgebra
include("../../../objects/mesh.jl")
using .IRMesh
using Firefly
cfg = Firefly.Config

# Config parameters
outdir = cfg.outdir
prefix = cfg.prefix
filetype = cfg.filetype
nx, ny, nz = cfg.k_mesh
dim = cfg.dimension
if dim == 2
    nz = 1
end
nk = nx * ny * nz
mu = cfg.fermi_energy
if cfg.mu_from_n
    n_field = Field_R(outdir*prefix*"_E_vs_n."*filetype)
    mu = n_field(cfg.num_electrons)
    println("Shifting mu to $(mu) based on electron number $(cfg.num_electrons)")
end
BZ = cfg.brillouin_zone
beta = 1 / cfg.Temperature
Z = cfg.qp_weight

function get_kvec(ix, iy, iz, nx, ny, nz)
    kvec = [ix / nx - 0.0, iy / ny - 0.0, iz / nz - 0.0]
    if dim < 3
        kvec[3] = 0.0
    elseif dim < 2
        kvec[2] = 0.0
    end
    kvec = BZ * kvec
    return kvec
end

function fill_energy_mesh(band)
    ek = Array{Float32}(undef, nx, ny, nz)
    for i in 1:nx, j in 1:ny, k in 1:nz
        kvec = get_kvec(i - 1, j - 1, k - 1, nx, ny, nz)
        ek[i, j, k] = band(kvec)
    end
    return ek
end

function fill_sigma_mesh(sigma, iw)
    nw = length(iw)
    Ekw = Array{ComplexF32}(undef, nw, nx, ny, nz)
    for i in 1:nx, j in 1:ny, k in 1:nz, l in 1:nw
        kvec = get_kvec(i - 1, j - 1, k - 1, nx, ny, nz)
        Ekw[i, j, k] = sigma(kvec, imag(iw[l]))
    end
    return Ekw
end

function main()
    println("Calculating Response Function χ(iν,k)")

    # Load band structure
    println("Constructing Bands")
    band = Bands()

    # Fill energy mesh
    println("Filling Energy Mesh")
    ek = fill_energy_mesh(band)
    minval, maxval = minimum(ek), maximum(ek)
    println("Energy range: [$minval, $maxval]")
    D = maxval - minval

    # Create IR mesh
    mesh = IR_Mesh(D)
    iw, iv = get_iw_iv(mesh)

    # Construct Green's function G(iω,k) = 1/(iω - ε(k) + μ)
    println("Constructing G(iω,k)")
    ek_reshaped = reshape(ek, 1, nx, ny, nz)
    iw_reshaped = reshape(iw, mesh.fnw, 1, 1, 1)

    #filename = outdir * prefix * "_self_energy.h5"
    #println("Checking if Self-Energy file `$filename` exists.")
    #Sigma = Field_C(filename)
    #Ekw = fill_sigma_mesh(Sigma, iw)
    
    #Gkw = 1.0 ./ (iw_reshaped .- (ek_reshaped .- mu))
    Gkw = Z ./ (iw_reshaped .- (Z .* ek_reshaped .- mu))

    # Calculate χ(iν,k) via convolution: χ = -G(iω,k) * G(-iω,-k)
    # Transform to (r,τ) space
    println("Transforming to real space and imaginary time")
    Grt = kw_to_rtau(Gkw, 'F', mesh)

    # Original convolution formula (testing with fixed Hamiltonian)
    println("Computing convolution χ(r,τ) = G(r,τ) × G(r,-τ)")
    chi_rt = Grt .* reverse(Grt, dims=1)

    # Transform back to (k,ν)
    println("Transforming back to momentum and bosonic frequency")
    chi_kw = rtau_to_kw(chi_rt, 'B', mesh)

    # Center k-points correctly: shift from (0,0) to (2π,2π) -> (-π,-π) to (π,π)
    println("Centering k-points")
    for i in 1:mesh.bnw
        chi_kw[i, :, :, :] .= fftshift(chi_kw[i, :, :, :])
    end

    println("Max χ: $(maximum(abs.(chi_kw)))")
    println("Min χ: $(minimum(abs.(chi_kw)))")

    # Save result
    kmesh = cfg.k_mesh
    BZ_in = BZ
    if dim == 2
        kmesh = kmesh[1:end-1]
        BZ_in = BZ_in[1:end-1, 1:end-1]
    end

    output_file = outdir * prefix * "_chi." * filetype
    println("Saving χ(iν,k) to $output_file")
    save_data!(output_file, chi_kw, kmesh, BZ_in, w_points = imag.(iv))

    return maximum(abs.(chi_kw))
end

function run()
    max_chi = main()
    return max_chi
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end
