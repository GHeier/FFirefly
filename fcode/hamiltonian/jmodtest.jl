const testmodule = "/home/g/Research/fcode/fcode/hamiltonian/jmodtest.so"

function epsilon(n::Int, k::Vector{Float64})
    return ccall((:epsilon_julia, testmodule), Float64, (Int, Ptr{Float64}, Int), n, k, length(k))
end

function load_config(path::String)
    ccall((:load_config_julia, testmodule), Cvoid, (Cstring,), path)
end

path = "/home/g/Research/fcode/build/bin/input.cfg"
load_config(path)
println(epsilon(1, [0.1, 0.2, 0.3]))
