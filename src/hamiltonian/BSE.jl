module BSE
include("../objects/mesh.jl")
using .IRMesh
using Firefly
cfg = Firefly.Config

nx, ny, nz = cfg.k_mesh
dim = cfg.dimension
if dim == 2
    nz = 1
end
nbnd = cfg.nbnd
mu = cfg.fermi_energy
U = cfg.onsite_U
BZ = cfg.brillouin_zone

function BSE_node()
    band = Bands()
    ek = fill_energy_mesh(band)
    mesh = IR_Mesh(ek)
    iw, iv = get_iw_iv(mesh)
    fnw, bnw, fntau, bntau = mesh.fnw, mesh.bnw, mesh.fntau, mesh.bntau
    iw = reshape(iw, fnw, 1, 1, 1)
    # Preallocate all arrays
    Sigma = zeros(ComplexF64, fnw, nx, ny, nz)
    G     = zeros(ComplexF64, fnw, nx, ny, nz)
    G_rt  = Array{ComplexF64}(undef, fntau, nx, ny, nz)
    V     = Array{ComplexF64}(undef, bnw, nx, ny, nz)
    V_rt  = Array{ComplexF64}(undef, bntau, nx, ny, nz)

    max_iters = 2
    for _ in 1:max_iters
        @. G = 1 / (iw - ek + mu - Sigma)

        kw_to_rtau!(G_rt, G, 'F', mesh)
        @. G_rt = G_rt * G_rt

        rtau_to_kw!(V, G_rt, 'B', mesh)

        @. V = (-3/2) * (U^2 * V / (1 - U * V)) + (1/2) * (U^2 * V / (1 + U * V))

        kw_to_rtau!(V_rt, V, 'B', mesh)
        @. V_rt = V_rt * G_rt

        rtau_to_kw!(Sigma, V_rt, 'F', mesh)
    end
end


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


function fill_energy_mesh(band)
    ek = Array{Float64}(undef, 1, nx, ny, nz)
    for i in 1:nx, j in 1:ny, k in 1:nz
        kvec = get_kvec(i, j, k, nx, ny, nz)
        ek[1, i, j, k] = band(1, kvec)
    end
    return ek
end


end # module
