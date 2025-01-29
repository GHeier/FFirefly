module Eliashberg

t1 = time()
include("../objects/fields.jl")
using .Fields
include("../objects/mesh.jl")
using .IRMesh

using CUDA
using Plots
using SparseIR
import SparseIR: Statistics, value, valueim
using Base.Threads, JLD
using LinearAlgebra, Printf, PyCall
np = pyimport("numpy")

t2 = time()
println("Time to load: ", t2 - t1)

fcode = pyimport("fcode")
t3 = time()
println("Time to load fcode: ", t3 - t2)
cfg = fcode.config

np_BZ = np.array(cfg.brillouin_zone)

nx, ny, nz = cfg.k_mesh
dim = cfg.dimension
if dim == 2
    nz = 1
end
nw = cfg.w_pts 
U = cfg.onsite_U
BZ = cfg.brillouin_zone
mu = cfg.fermi_energy
t4 = time()
println("Time to load config: ", t4 - t3)

#const beta = 1/cfg.Temperature
const beta = 1000.0
println("Beta: ", beta)
const pi = π



function to_IBZ(k)
    tolerance = 1e-10  # small tolerance to account for floating-point errors
    q = abs.(k)
    q .= ifelse.(abs.(q .- 2π) .< tolerance, 0.0, ifelse.(q .> π, -(q .- 2π), q))
    q = abs.(q)
    q .+= 1e-4
    q .= ifelse.(q .> 2π, q .- 1e-4, q)
    return q
end

function epsilon(k)
    if dim == 2
        return -2 * (cos(k[1]) + cos(k[2]))
    end
    return -2 * (cos(k[1]) + cos(k[2]) + cos(k[3]))
end

function V(k1, k2, w)
    qm = k1 - k2
    qm = to_IBZ(qm)
    qm .= round.(qm, digits=6)
    iw = round(imag(w), digits=6)
    Xm = Susceptibility(qm[1], qm[2], qm[3], iw)
    Vm = U^2 * Xm / (1 - U*Xm) + U^3 * Xm^2 / (1 - U^2 * Xm^2)
    return Vm
end

function get_bandwidth()
    maxval = -1000
    minval = 1000
    for i in 1:200, j in 1:200, k in 1:200
        kvec = BZ * [i / 200, j / 200, k / 200]
        eps = epsilon(kvec)
        maxval = max(maxval, eps)
        minval = min(minval, eps)
    end
    return maxval - minval
end

function paper_V(w)
    lambda = 2.0
    return lambda * 0.01 / (imag(w)^2 + 0.01)
end

function paper2_V(w)
    n = (beta * imag(w) / pi - 1) / 2
    lambda = 2.0
    v = beta * 1.5 / (2 * pi)
    return -2 * lambda * v^2 / (n^2 + v^2)
end

function fill_V_arr!(V_arr, iw)
    Vw_arr = Array{Complex{Float64}}(undef, bnw, nx, ny)
    for i in 1:nx, j in 1:ny, k in 1:nz, l in 1:bnw
        kvec = BZ * [i / nx, j / ny, k / nz]
        w1 = iw[l]
        #Vw_arr[l, i, j] = V(kvec, [0,0,0], w1)
        Vw_arr[l, i, j] = paper2_V(w1)
    end
    temp = k_to_r(mesh, Vw_arr)
    V_arr .= wn_to_tau(mesh, Bosonic(), temp)
    println("Filled V_arr")
end

function initialize_phi_Z_chi!(phi_arr, Z_arr, chi_arr, iw)
    for i in 1:fnw, j in 1:nx, k in 1:ny
        w = iw[i]
        kvec = BZ * [j / nx, k / ny, 0 / nz]
        dwave = (cos(kvec[1]) - cos(kvec[2])) / 2
        phi_arr[i, j, k] = 1 / abs(imag(w)^2 + 1)# * dwave
        Z_arr[i, j, k] = 1.0 
        chi_arr[i, j, k] = 0.0
    end
    println("Initialized phi and Z")
end

function condense_to_F!(phi, Z, F_arr, Fz_arr, iw)
    for i in 1:nx, j in 1:ny, k in 1:nz, l in 1:fnw
        w = iw[l]
        kvec = BZ * [i / nx, j / ny, k / nz]
        denom = get_denominator(w, phi[l, i, j], Z[l, i, j], 2 - mu)
        F_arr[l, i, j] = phi[l, i, j] / denom
        Fz_arr[l, i, j] = Z[l, i, j] * w / denom / 1.0im
    end
end

function condense_to_F_and_G!(phi, Z, chi, F_arr, G_arr, iw)
    for i in 1:nx, j in 1:ny, k in 1:nz, l in 1:fnw
        w = iw[l]
        kvec = BZ * [i / nx, j / ny, k / nz]
        denom = get_denominator(w, phi[l, i, j], Z[l, i, j], chi[l, i, j] + epsilon(kvec) - mu)
        F_arr[l, i, j] = -phi[l, i, j] / denom
        G_arr[l, i, j] = -(w * Z[l, i, j] + epsilon(kvec) - mu + chi[l, i, j]) / denom
    end
end

function update!(F, G, V_arr, phi, Z, chi, iw, sigma)
    temp = k_to_r(mesh, F)
    F_rt = wn_to_tau(mesh, Fermionic(), temp)
    temp = k_to_r(mesh, G)
    G_rt = wn_to_tau(mesh, Fermionic(), temp)

    phit = V_arr .* F_rt
    sigmat = -V_arr .* G_rt

    temp = r_to_k(mesh, phit)
    phi .= tau_to_wn(mesh, Fermionic(), temp)
    temp = r_to_k(mesh, sigmat)
    sigma .= tau_to_wn(mesh, Fermionic(), temp)

    for i in 1:fnw, j in 1:nx, k in 1:ny
        w = iw[i]
        sp, sm = sigma[i, j, k], sigma[fnw - i + 1, j, k]
        Z[i, j, k] = 1.0 - 0.5 * (sp - sm) / w
        chi[i, j, k] = 0.5 * (sp + sm)
    end
end

function update_phi_Z!(F_arr, Fz_arr, V_arr, phi_arr, Z_arr, iw)
    #temp = k_to_r(mesh, F_arr)
    F_rt = wn_to_tau(mesh, Fermionic(), F_arr)
    #temp = k_to_r(mesh, Fz_arr)
    Fz_rt = wn_to_tau(mesh, Fermionic(), Fz_arr)

    phit = V_arr .* F_rt
    Zt = V_arr .* Fz_rt

    #temp = r_to_k(mesh, phit)
    phi_arr .= tau_to_wn(mesh, Fermionic(), phit)
    #temp = r_to_k(mesh, Zt)
    Z_arr .= tau_to_wn(mesh, Fermionic(), Zt)

    for i in 1:fnw
        w = iw[i]
        Z_arr[i, :, :] .= 1 .- 1.0 / imag(w) .* Z_arr[i, :, :]
    end
    phi_arr .= phi_arr .* (-1.0 )
end

function test_update_phi_Z(F_arr, Fz_arr, iw)
    phi_new = Array{Complex{Float64}}(undef, fnw, nx, ny)
    Z_new = Array{Complex{Float64}}(undef, fnw, nx, ny)
    @threads for i in 1:fnw
        for j in 1:nx, k in 1:ny
            sum_phi = 0.0 + 0.0im
            sum_Z = 0.0 + 0.0im
            k1 = BZ * [j / nx, k / ny, 0 / nz]
            w1 = iw[i]
            for l in 1:fnw, m in 1:nx, n in 1:ny
                k2 = BZ * [m / nx, n / ny, 0 / nz]
                w2 = iw[l]
                F = F_arr[l, m, n]
                Fz = Fz_arr[l, m, n]
                #V_val = V(k1, k2, w1)
                V_val = paper_V(w1 - w2)
                sum_phi += V_val * F
                sum_Z += V_val * Fz
            end
            phi_new[i, j, k] = sum_phi / beta
            Z_new[i, j, k] = 1 - sum_Z / (beta * w1)
        end
    end
    return phi_new, Z_new
end

function newfunc()
    println("Beginning Eliashberg")
    wmax = get_bandwidth()
    println("Bandwidth: ", wmax)
    IR_tol = 1e-15
    scf_tol = 1e-4
    IR_basis_set = FiniteTempBasisSet(beta, Float64(wmax), IR_tol)
    global mesh = Mesh(IR_basis_set, nx, ny)
    global fnw = mesh.fnw
    global fntau = mesh.fntau
    global bnw = mesh.bnw
    global bntau = mesh.bntau
    println(fnw, " ", fntau, " ", bnw, " ", bntau)
    println("IRMesh created")
    iw = Array{ComplexF64}(undef, fnw)
    iv = Array{ComplexF64}(undef, bnw)
    npts = Array{Integer}(undef, fnw)
    for i in 1:fnw iw[i] = valueim(mesh.IR_basis_set.smpl_wn_f.sampling_points[i], beta) end
    for i in 1:bnw iv[i] = valueim(mesh.IR_basis_set.smpl_wn_b.sampling_points[i], beta) end
    for i in 1:fnw npts[i] = round((iw[i].im * beta / pi - 1) / 2) end 

    println("Filled iw and iv")
    println("Min & Max iw: ", minimum(imag.(iw)), " ", maximum(imag.(iw)))
    input_data_file = "/home/g/Research/fcode/chi_mesh_dynamic.dat"
    #input_data_file = "/home/g/Research/fcode/ckio_ir.dat"
    global Susceptibility = Fields.create_interpolation(input_data_file, dims=4, field_type=:scalar, value_type=ComplexF64)

    phi_arr = Array{ComplexF64}(undef, fnw, nx, ny)
    Z_arr = Array{ComplexF64}(undef, fnw, nx, ny)
    chi_arr = Array{ComplexF64}(undef, fnw, nx, ny)
    println("Initializing phi, Z, and chi")
    initialize_phi_Z_chi!(phi_arr, Z_arr, chi_arr, iw)

    println("maxphi: ", maximum(abs.(phi_arr)))
    println("maxZ: ", maximum(abs.(Z_arr)))
    println("maxchi: ", maximum(abs.(chi_arr)))

    sigma = Array{ComplexF64}(undef, fnw, nx, ny)
    V_arr = Array{ComplexF64}(undef, bntau, nx, ny)
    fill_V_arr!(V_arr, iv)
    println("maxV: ", V_arr[argmax(abs.(V_arr))])
    println("maxV: ", maximum(abs.(V_arr)))

    F_arr = Array{ComplexF64}(undef, fnw, nx, ny)
    G_arr = Array{ComplexF64}(undef, fnw, nx, ny)
    #Fz_arr = Array{ComplexF64}(undef, fnw, nx, ny)
    #condense_to_F!(phi_arr, Z_arr, F_arr, Fz_arr, iw)
    condense_to_F_and_G!(phi_arr, Z_arr, chi_arr, F_arr, G_arr, iw)
    println("maxF: ", F_arr[argmax(abs.(F_arr))])
    println("maxG: ", G_arr[argmax(abs.(G_arr))])
    #println("maxFz: ", Fz_arr[argmax(abs.(Fz_arr))])


    phierr = 0.0
    prev_phi_arr = copy(phi_arr)
    #phi_test, Z_test = Array{ComplexF64}(undef, fnw, nx, ny), Array{ComplexF64}(undef, fnw, nx, ny)
    println("Starting Convergence Loop")
    iterations = 300
    for i in 1:iterations
        update!(F_arr, G_arr, V_arr, phi_arr, Z_arr, chi_arr, iw, sigma)
        #F_test = copy(F_arr)
        #Fz_test = copy(Fz_arr)
        #update_phi_Z!(F_arr, Fz_arr, V_arr, phi_arr, Z_arr, iw)
        #phi_test, Z_test = test_update_phi_Z(F_test, Fz_test, iw)

        phierr = round(sum(abs.(prev_phi_arr - phi_arr)) / (fnw * nx * ny),digits=8)
        print("Iteration $i: Error = $phierr           \r")

        if abs(phierr) < scf_tol || maximum(abs.(phi_arr)) < 1e-10
            println("Converged after ", i, " iterations              ")
            break
        end
        prev_phi_arr = copy(phi_arr)
        condense_to_F_and_G!(phi_arr, Z_arr, chi_arr, F_arr, G_arr, iw)
        #condense_to_F!(phi_arr, Z_arr, F_arr, Fz_arr, iw)
        #condense_to_F!(phi_test, Z_test, F_test, Fz_test, iw)
    end

    #phi_w, Z_w = test_update_phi_Z(F_arr, Fz_arr, iw)
    max_phi = maximum(abs.(phi_arr))
    max_Z = maximum(abs.(Z_arr))

    println("Error: ", phierr)
    println("Max phi: ", max_phi)
    println("Max Z: ", max_Z)
    #gap = phi_arr ./ Z_arr
    #gap_w = phi_FS_average(gap)
    #p = plot(imag.(iw), gap_w)
    #display(p)
    #readline()
    phi_w = phi_FS_average(phi_arr)
    p = plot(npts, phi_w, xlimits=(-500, 500))
    display(p)
    readline()
    Z_w = phi_FS_average(Z_arr)
    p = plot(npts, Z_w, xlimits=(-500, 500))
    display(p)
    readline()
    sigma_w = phi_FS_average(sigma)
    p = plot(npts, sigma_w, xlimits=(-500, 500))
    display(p)
    readline()
    p = plot(imag.(iw), ylimits=(-5, 5))
    display(p)
    readline()
    println(imag.(iw))
end 

function phi_FS_average(phi)
    phi_avg = Array{Float64}(undef, fnw)
    for i in 1:fnw
        phi_avg[i] = 0.0 + 0.0im
        for j in 1:nx, k in 1:ny
            phi_avg[i] += phi[i, j, k].re^2
        end
        phi_avg[i] /= nx * ny
    end
    return phi_avg
end

function get_denominator(w1, phi_el, Z_el, eps)
    return abs2(w1 * Z_el) + abs2(eps) + abs2(phi_el)
end


function k_integral(k, w, w1, phi, Z)
    phi_int = Z_int = 0.0 + 0.0im
    n = trunc(Int, (imag(w) * beta / pi - 1) / 2)
    temp = 0

    for i in 1:nx 
        for j in 1:ny 
            for l in 1:nz
                k1 = BZ * [i / nx, j / ny, l / nz]
                phi_el, Z_el= phi[i, j, l, n], Z[i, j, l, n]
                V_val = (V(k, k1, w) + V(k, -k1, w)) / 2.0
                denom = get_denominator(w1, phi_el, Z_el, epsilon(k1) - mu)
                temp = max(temp, denom)
                phi_int += V_val * phi_el / denom
                Z_int += V_val * w1 / w * Z_el / denom
            end
        end
    end

    #println(w.im, " ", w1.im, " ", Z_int.re, " ", temp)
    return phi_int / (nx * ny * nz), Z_int / (nx * ny * nz)
end

function eliashberg_sum(phi, Z)
    new_phi = zeros(Complex{Float64}, size(phi))
    new_Z = zeros(Complex{Float64}, size(Z))

    for a in 1:nx, b in 1:ny, c in 1:nz
        print("Iteration ", (a-1)*nx^(dim-1) + (b-1)*ny^(dim-2) + c-1, " out of ", nx*ny*nz, "      \r")
        flush(stdout)
        @threads for i in 1:nw 
            phi_sum = Z_sum = 0.0 + 0.0im
            for j in 1:nw
                w = Complex(0.0, pi*(2*(i-1) + 1) / beta)
                w1 = Complex(0.0, pi*(2*(j-1) + 1) / beta)
                k = BZ * [a / nx, b / ny, c / nz]
                phi_int, Z_int = k_integral(k, w, w1, phi, Z)
                phi_sum += phi_int
                Z_sum += Z_int
            end
            new_phi[a, b, c, i] -= phi_sum / beta
            new_Z[a, b, c, i] += (1 - Z_sum / beta)
            #println("phi_sum: ", phi_sum, " Z_sum: ", Z_sum)
        end
    end
    return new_phi, new_Z
end

function evaluate_eliashberg()
    phi = Z = ones(Complex{Float64}, nx, ny, nz, nw)
    iterations = 100

    input_data_file = "/home/g/Research/fcode/chi_mesh_dynamic.dat"
    global Susceptibility = Fields.create_interpolation(input_data_file, dims=4, field_type=:scalar, value_type=ComplexF64)
    value = Susceptibility(0.005, 0.005, 0.005, 0.005)
    println("Value: ", value)
    #fcode.load_global_chi(input_data_file)

    println("Starting Eliashberg calculation")
    prev_phi = prev_Z = 0.0
    for i in 0:iterations
        phi, Z = eliashberg_sum(phi, Z)
        max_phi = maximum(abs.(phi))
        max_Z = maximum(abs.(Z))
        println("\nMax phi: ", max_phi)
        println("Max Z: ", max_Z)
        error = abs(prev_phi - max_phi) + abs(prev_Z - max_Z)
        println("Error: ", maximum(error))
        if maximum(error) < 1e-4
            println("Converged after ", i, " iterations")
            break
        end
        prev_phi, prev_Z = max_phi, max_Z
    end
    save("phi.jld", "phi", phi)
    save("Z.jld", "Z", Z)
end

function multiply_slow(V, f)
    result = Array{ComplexF64}(undef, size(f))
    fnw = length(f)
    for i in 1:fnw
        result[i] = 0.0 + 0.0im
        for j in 1:fnw
            result[i] += V[i, j] * f[j]
        end
    end
    return result / beta
end

function multiply(f, V, mesh)
    fntau = mesh.fntau        
    #fntau = 100

    ft = wn_to_tau(mesh, Fermionic(), f)
    Vt = wn_to_tau(mesh, Bosonic(), V)
    #revF = Array{ComplexF64}(undef, fntau)
    #for i in 1:fntau
    #    revF[i] = f[fntau - i + 1]
    #end
    resultt = Vt .* ft
    result = tau_to_wn(mesh, Fermionic(), resultt)

    #temp = k_to_r(mesh, f)
    #t_result = temp .* temp
    #result = tau_to_wn(mesh, r_result)

    return result
end

function test_iw()
    IR_basis_set = FiniteTempBasisSet(beta, Float64(10), 1e-10)
    basis = Mesh(IR_basis_set, 10)
    fnw = basis.fnw
    fntau = basis.fntau
    bnw = basis.bnw
    bntau = basis.bntau
    iw = Array{ComplexF64}(undef, fnw)
    iv = Array{ComplexF64}(undef, bnw)
    for i in 1:fnw
        iw[i] = valueim(basis.IR_basis_set.smpl_wn_f.sampling_points[i], beta)
        if i <= bnw iv[i] = valueim(basis.IR_basis_set.smpl_wn_b.sampling_points[i], beta) end
    end
    println("fnw: ", fnw)
    func(x) = 1 / (abs(x) - 1)
    f = Array{ComplexF64}(undef, fnw)
    V = Array{ComplexF64}(undef, bnw)
    for i in 1:fnw
        w = iw[i]
        f[i] = func(w)
        if i<=bnw V[i] = paper_V(iv[i]) end
    end
    V_arr = Array{ComplexF64}(undef, bnw, bnw)
    for i in 1:bnw, j in 1:bnw
        w1 = iv[i]
        w2 = iv[j]
        V_arr[i, j] = paper_V(w1 - w2)
    end
    result = multiply(f, V, basis)
    slow_result = multiply_slow(V_arr, f)

    p = plot(imag.(iw), real.(result))
    plot!(imag.(iw), imag.(result))
    plot!(imag.(iw), real.(slow_result))
    plot!(imag.(iw), imag.(slow_result))
    display(p)
    readline()
end

function f1(x)
    return exp(-x^2)
    return 1 / (sin(x) + 1)
    return sin(x)
end

function test_multiply()
    pts = 1000
    IR_basis_set = FiniteTempBasisSet(beta, Float64(20), 1e-10)
    amesh = Mesh(IR_basis_set, pts)
    bnw = amesh.bnw
    bntau = amesh.bntau
    wpts = 100
    println("iw size: ", wpts)
    Vw = Array{ComplexF64}(undef, wpts, pts)
    phiw = Array{ComplexF64}(undef, wpts, pts)
    for i in 1:pts
        x = 2π * i / (pts-1)
        for j in 1:wpts
            w = (pi*(2*(j-1) + 1) / beta)^2im
            Vw[j, i] = f1(x)
            phiw[j, i] = f1(x)
        end
    end
    V = k_to_r(amesh, Vw)
    #V = wn_to_tau(amesh, Bosonic(), temp)
    phi = k_to_r(amesh, phiw)
    #phi = wn_to_tau(amesh, Bosonic(), temp)

    println("size of V: ", size(V))
    println("size of phi: ", size(phi))
    flipV = Array{ComplexF64}(undef, wpts, pts)
    for i in 1:pts, j in 1:wpts
        flipV[j, i] = V[wpts + 1 - j, i]
    end

    result = V .* phi

    newphi = r_to_k(amesh, result)
    #newphi = tau_to_wn(amesh, Bosonic(), temp)
    println("Completed multiplication")
    println(newphi[1, 1])
    p = plot(0:1/(pts-1):1, real.(newphi[1,:]))
    plot!(0:1/(pts-1):1, imag.(newphi[1,:]))
    display(p)
    readline()
    #p = plot(0:1/(wpts-1):1, real.(newphi[:,1]))
    #plot!(0:1/(wpts-1):1, imag.(newphi[:,1]))
    #display(p)
    #readline()
end

function test_multiply2()
    pts = 1000
    #IR_basis_set = FiniteTempBasisSet(beta, Float64(20), 1e-10)
    #amesh = Mesh(IR_basis_set, pts)
    #bnw = amesh.bnw
    wpts = 1
    println("iw size: ", wpts)
    W = Array{ComplexF64}(undef, wpts, pts)
    f = Array{ComplexF64}(undef, wpts, pts, wpts, pts)
    g = Array{ComplexF64}(undef, wpts, pts)
    for i in 1:pts, j in 1:pts, k in 1:wpts, l in 1:wpts
        x = 2π * i / (pts-1)
        y = 2π * j / (pts-1)
        w1 = (pi*(2*(j-1) + 1) / beta)^2im
        w2 = (pi*(2*(l-1) + 1) / beta)^2im
        f[k, i, l, j] = f1(x - y)
    end
    println("Filled f")
    for i in 1:pts, j in 1:wpts
        x = 2π * i / (pts-1)
        w = (pi*(2*(j-1) + 1) / beta)^2im
        g[j, i] = f1(x)
    end
    for i in 1:pts, j in 1:wpts
        sum = 0.0
        for k in 1:pts, l in 1:wpts
            sum += f[j, i, l, k] * g[l, k]
        end
        W[j, i] = sum / (pts*wpts)
    end
    p = plot(0:1/(pts-1):1, real.(W[1,:]))
    plot!(0:1/(pts-1):1, imag.(W[1,:]))
    display(p)
    println(W[1,1])
    readline()
    #p = plot(0:1/(wpts-1):1, real.(W[:,1]))
    #plot!(0:1/(wpts-1):1, imag.(W[:,1]))
    #display(p)
    #readline()
end

function restart()
    println("Restarting Eliashberg")
    println(Threads.nthreads(), " threads available")
    IR_tol = 1e-10
    scf_tol = 1e-4
    IR_basis_set = FiniteTempBasisSet(beta, Float64(8), IR_tol)
    basis = Mesh(IR_basis_set, nx, ny)
    fnw = basis.fnw
    fnw = 1000

    phi = Array{ComplexF64}(undef, fnw)
    Z = Array{ComplexF64}(undef, fnw)
    npts = Array{Integer}(undef, fnw)
    iw = Array{ComplexF64}(undef, fnw)

    for i in 1:fnw
        npts[i] = -fnw/2 + i - 1
        #npts[i] = i - 1
        #iw[i] = valueim(basis.IR_basis_set.smpl_wn_f.sampling_points[i], beta) 
        iw[i] = Complex(0.0, (pi*(2*(npts[i] - 1) + 1) / beta))
        phi[i] = 1 / abs(iw[i] + 1)
        Z[i] = 1 / abs(iw[i] + 1)
    end

    #V(x) = 2 * (beta * 1.5 / (2*pi))^2 / (x^2 + (beta * 1.5 / (2*pi))^2)
    w0 = 1e-1
    #V(x) = abs(x) < w0 ? -1.0 : 0.0
    V(x) = paper2_V(x)

    for iters in 1:100
        print("Iteration $iters \r")
        new_phi = zeros(Complex{Float64}, fnw)
        new_Z = ones(Complex{Float64}, fnw)
        @threads for i in 1:fnw
            for j in 1:fnw, k in 1:nx, l in 1:ny
                kvec = BZ * [k / nx, l / ny, 0 / nz]
                denom = get_denominator(iw[j], phi[j], Z[j], epsilon(kvec) - mu)
                n = npts[i] - npts[j]
                #nn = npts[i] + npts[j] + 1
                new_phi[i] -= 1 / beta * V(n)* phi[j] / denom
                new_Z[i] -= 1 / beta * V(n) * (iw[j] / iw[i]) * Z[j] / denom
            end
        end
        if maximum(abs.(new_phi - phi)) < scf_tol
            println("Converged after ", iters, " iterations")
            break
        end
        phi = new_phi
        Z = new_Z
    end
    println("Max phi: ", maximum(abs.(phi)))
    println("Max Z: ", maximum(abs.(Z)))
    p = plot(npts, real.(phi))
    display(p)
    readline()
    p = plot(npts, real.(Z))
    display(p)
    readline()
end

function gpu_sum()
    # Initialize data on GPU
    scf_tol = 1e-4
    fnw = 1000
    npts = Array{Integer}(undef, fnw)
    iw = Array{ComplexF64}(undef, fnw)
    phi = Array{ComplexF64}(undef, fnw)
    Z = Array{ComplexF64}(undef, fnw)
    for i in 1:fnw
        npts[i] = -fnw/2 + i - 1
        iw[i] = Complex(0.0, (pi*(2*(npts[i] - 1) + 1) / beta))
        phi[i] = 1 / abs(iw[i] + 1)
        Z[i] = 1 / abs(iw[i] + 1)
    end

    iw_gpu = CuArray(iw)
    phi_gpu = CuArray(phi)
    Z_gpu = CuArray(Z)

    BZ_gpu = CuArray(BZ)

    for iters in 1:20
        print("Iteration $iters: ")
        new_phi = CuArray(zeros(Complex{Float64}, fnw))
        new_Z = CuArray(ones(Complex{Float64}, fnw))

        #@device_code_warntype kernel_summation!(
        #    new_phi, new_Z, phi_gpu, Z_gpu, iw_gpu, nx, ny, beta, BZ_gpu, mu
        #)
        # GPU Kernel for summation
        @cuda threads=256 blocks=cld(fnw, 256) kernel_summation!(
             new_phi, new_Z, phi_gpu, Z_gpu, iw_gpu, nx, ny, beta, BZ_gpu, mu
        )

        # Convergence check
        error = CUDA.maximum(abs.(new_phi - phi_gpu))
        println("Error = ", error)
        if error < scf_tol
            println("Converged after $iters iterations")
            break
        end

        # Update variables
        phi_gpu .= new_phi
        Z_gpu .= new_Z
    end
    phi = Array(phi_gpu)
    Z = Array(Z_gpu)
    println("Max phi: ", maximum(abs.(phi)))
    println("Max Z: ", maximum(abs.(Z)))
    p = plot(npts, real.(phi))
    display(p)
    readline()
    p = plot(npts, real.(Z))
    display(p)
    readline()
end

# GPU Kernel Definition
function kernel_summation!(
    new_phi::CuDeviceVector{ComplexF64},
    new_Z::CuDeviceVector{ComplexF64},
    phi::CuDeviceVector{ComplexF64},
    Z::CuDeviceVector{ComplexF64},
    iw::CuDeviceVector{ComplexF64},
    nx::Int,
    ny::Int,
    beta::Float64,
    BZ::CuDeviceMatrix{Float64},
    mu::Float64
)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    if i > length(new_phi)
        return
    end

    temp_phi = ComplexF64(0.0, 0.0)
    temp_Z = ComplexF64(0.0, 0.0)

    V(x) = 1.0 / (imag(x)^2 + 0.01)
    e(kx, ky, kz) = -2 * (cos(kx) + cos(ky))
    get_denom(iw, phi, Z, eps) = abs2(imag(iw)*Z) + abs2(eps) + abs2(phi)

    for j in 1:length(phi)
        for k in 1:nx, l in 1:ny
            kx = BZ[1, 1] * (k / nx) + BZ[1, 2] * (l / ny) + BZ[1, 3] * 0.0
            ky = BZ[2, 1] * (k / nx) + BZ[2, 2] * (l / ny) + BZ[2, 3] * 0.0
            kz = BZ[3, 1] * (k / nx) + BZ[3, 2] * (l / ny) + BZ[3, 3] * 0.0
            denom = get_denom(iw[j], phi[j], Z[j], e(kx, ky, kz) - mu)
            w = iw[i] - iw[j]

            V_val = ComplexF64(-1.0, 0.0)
            V_val = V(w)
            # Accumulate results
            temp_phi -= (1 / beta) * V_val * phi[j] / denom
            temp_Z -= (1 / beta) * V_val * (iw[j] / iw[i]) * Z[j] / denom
        end
    end

    new_phi[i] = temp_phi
    new_Z[i] = temp_Z
    return nothing
end

end # module


#using .Eliashberg
#println(Threads.nthreads(), " threads available")
#Eliashberg.evaluate_eliashberg()
Eliashberg.newfunc()
#Eliashberg.restart()
#Eliashberg.gpu_sum()

#Eliashberg.test_iw()
#Eliashberg.test_multiply()
#Eliashberg.test_multiply2()

