#  g++ -shared -fPIC ../../build/CMakeFiles/fcode.x.dir/fcode/config/load/c_config.c.o ../../build/CMakeFiles/fcode.x.dir/fcode/objects/vec.cpp.o ../../build/CMakeFiles/fcode.x.dir/fcode/hamiltonian/band_structure.cpp.o ../../build/CMakeFiles/fcode.x.dir/fcode/config/load/cpp_config.cpp.o jmodtest.o -o jmodtest.so
const testmod = "/home/g/Research/fcode/fcode/hamiltonian/jmodtest.so"

function epsilon(n::Int, k::Vector{Float64})
    return ccall((:epsilon_julia, testmod), Float64, (Int, Ptr{Float64}, Int), n, k, length(k))
end

function load_config!(path::String)
    ccall((:load_config_julia, testmod), Cvoid, (Cstring,), path)
end

# Create a new CMF instance
function create_CMF_CS()
    return ccall((:create_CMF_CS, testmod), Ptr{Cvoid}, ())
end

function cmf_cs_call(cmf::Ptr{Cvoid}, k::Vector{Float64}, w::Float64)
    result = Ref(ComplexF64(0.0, 0.0))  # Pass by reference
    ccall((:cmf_cs_call, testmod), ComplexF64, (Ptr{Cvoid}, Ptr{Float64}, Float64, Int, Ptr{ComplexF64}), cmf, k, w, length(k), result)
    return result[]
end

function cmf_cs_call2(cmf::Ptr{Cvoid}, w::Float64)
    result = Ref(ComplexF64(0.0, 0.0))  # Pass by reference
    ccall((:cmf_cs_call2, testmod), ComplexF64, (Ptr{Cvoid}, Float64, Ptr{ComplexF64}), cmf, w, result)
    return result[]
end

function load_cmf_cs(path::String)
    return ccall((:load_cmf_cs, testmod), Ptr{Cvoid}, (Cstring,), path)
end

function save_cmf_to_file!(path::String, cmf::Ptr{Cvoid})
    ccall((:cmf_save, testmod), Cvoid, (Cstring, Ptr{Cvoid}), path, cmf)
end

cmf = create_CMF_CS()
file = "/home/g/Research/Materials/Test/test_ckio_ir.dat"
cmf2 = load_cmf_cs(file)
println(cmf_cs_call(cmf2, [-3.14, -3.14, -3.14], 0.0))

#path = "/home/g/Research/fcode/build/bin/input.cfg"
#load_config!(path)
#println(epsilon(1, [0.1, 0.2, 0.3]))
