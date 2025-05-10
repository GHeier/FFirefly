module Test

using Base.Threads
using ThreadsX

using Firefly

cfg = Firefly.Config

nx, ny, nz = cfg.k_mesh
dim = cfg.dimension
if dim < 3
    nz = 1
end
BZ = cfg.brillouin_zone

outdir = cfg.outdir
prefix = cfg.prefix
wc = cfg.bcs_cutoff_frequency
mu = cfg.fermi_energy


function phi_energy_integral(N, phi, pts)
    de = 2 * wc / pts
    sum = ThreadsX.sum(1:pts) do i
        e = -wc + 2 * (i - 1) * wc / (pts - 1)
        N(e + mu) / (2 * (phi^2 + e^2)^(0.5)) 
    end
    return sum * de
end


function phi_energy_integral_convergence()
    N = Firefly.Field_R(outdir * prefix * "_DOS.dat")
    phi = 1.0
    pts = 100000
    iters = 200

    initial_integral = N(mu) * asinh(wc)
    final_phi = 2 * wc * exp(-1 / N(mu))

    #println("Expected initial_integral: ", initial_integral)
    for i in 1:iters
        integral_value = phi_energy_integral(N, phi, pts)
        dphi = -(1 - integral_value) * phi
        phi += dphi
        if (i - 1) % 40 == 0
            #println("Integral: ", integral_value, " phi: ", phi)
        end
        if abs(dphi / phi) < 1e-4
            #println("Integral: ", integral_value, " phi: ", phi)
            break
        end
    end
    println(phi)
    #println("Expected final integral: 1 expected final phi: ", final_phi)
end


function get_kvec(i, j, k)
    temp = BZ * [i / nx - 0.5, j / ny - 0.5, k / nz - 0.5]
    return temp[1:dim]
end


function get_kpts()
    kpts_list = Vector{Vector{Float64}}()
    for i in 1:nx, j in 1:ny, k in 1:nz
        kvec = get_kvec(i, j, k)
        e = epsilon(1, kvec)
        if abs(e - mu) > wc
            continue
        end
        push!(kpts_list, kvec)
    end
    return kpts_list
end


function phi_momentum_integral(phi, kpts)
    pts = length(kpts)
    #dk = 1 / pts
    sum = ThreadsX.sum(1:pts) do i
        k = kpts[i]
        e = epsilon(1, k) - mu
        if abs(e) < wc
            1 / (2 * (phi^2 + e^2)^(0.5)) 
        end
    end
    return sum / (nx * ny * nz)
end


function phi_momentum_integral_convergence()
    N = Firefly.Field_R(outdir * prefix * "_DOS.dat")
    phi = 1.0
    pts = 1000
    iters = 120

    initial_integral = N(mu) * asinh(wc)
    final_phi = 2 * wc * exp(-1 / N(mu))

    kpts = get_kpts()
    #println("Number of kpts: ", length(kpts))
    #println("Expected initial_integral: ", initial_integral)
    for i in 1:iters
        integral_value = phi_momentum_integral(phi, kpts)
        dphi = -(1 - integral_value) * phi
        phi += dphi
        if (i - 1) % 20 == 0
            #println("Integral: ", integral_value, " phi: ", phi)
        end
        if abs(dphi / phi) < 1e-4
            #println("Integral: ", integral_value, " phi: ", phi)
            break
        end
    end
    println(phi)
    #println("Expected final integral: 1 Expected final phi: ", final_phi)
end


function test()
    println("Num Threads: ", Threads.nthreads())
    phi_energy_integral_convergence()
    phi_momentum_integral_convergence()
end

end # module
