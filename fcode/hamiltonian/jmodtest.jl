#  g++ -shared -fPIC ../../build/CMakeFiles/fcode.x.dir/fcode/config/load/c_config.c.o ../../build/CMakeFiles/fcode.x.dir/fcode/objects/vec.cpp.o ../../build/CMakeFiles/fcode.x.dir/fcode/hamiltonian/band_structure.cpp.o ../../build/CMakeFiles/fcode.x.dir/fcode/config/load/cpp_config.cpp.o jmodtest.o -o jmodtest.so
const testmodule = "/home/g/Research/fcode/fcode/hamiltonian/jmodtest.so"

function epsilon(n::Int, k::Vector{Float64})
    return ccall((:epsilon_julia, testmodule), Float64, (Int, Ptr{Float64}, Int), n, k, length(k))
end

function load_config!(path::String)
    ccall((:load_config_julia, testmodule), Cvoid, (Cstring,), path)
end

# Create a new CMF instance
function create_CMF_CS()
    return ccall((:create_CMF_CS, testmodule), Ptr{Cvoid}, ())
end

function load_cmf_cs(path::String)
    return ccall((:load_cmf_cs, testmodule), Ptr{Cvoid}, (Cstring,), path)
end

function save_cmf_to_file!(path::String, cmf::Ptr{Cvoid})
    ccall((:cmf_save, testmodule), Cvoid, (Cstring, Ptr{Cvoid}), path, cmf)
end

cmf = create_CMF_CS()
file = "/home/g/Research/Materials/Test/test_ckio_ir.dat"
cmf2 = load_cmf_cs(file)

path = "/home/g/Research/fcode/build/bin/input.cfg"
load_config!(path)
println(epsilon(1, [0.1, 0.2, 0.3]))
