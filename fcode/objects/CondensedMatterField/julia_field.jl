#  g++ -shared -fPIC ../../build/CMakeFiles/fcode.x.dir/fcode/config/load/c_config.c.o ../../build/CMakeFiles/fcode.x.dir/fcode/objects/vec.cpp.o ../../build/CMakeFiles/fcode.x.dir/fcode/hamiltonian/band_structure.cpp.o ../../build/CMakeFiles/fcode.x.dir/fcode/config/load/cpp_config.cpp.o jmodtest.o -o jmodtest.so

#function epsilon(n::Int, k::Vector{Float64})
#    return ccall((:epsilon_julia, testmod), Float64, (Int, Ptr{Float64}, Int), n, k, length(k))
#end

module Field

#const testmod = "/home/g/Research/fcode/fcode/hamiltonian/jmodtest.so"
const testmod = "/home/g/Research/fcode/build/lib/libfcode.so"
mutable struct Field_C
    cmf::Ptr{Cvoid}
    function Field_C(cmf::Ptr{Cvoid})
        new(cmf)
    end

    function Field_C(filename::String)
        cmf = ccall((:load_cmf_cs, testmod), Ptr{Cvoid}, (Cstring,), filename)
        if cmf == C_NULL
            error("Failed to load CMF_CS from file: ", filename)
        end
        new(cmf)
    end

    function destroy!(csf::Field_C)
        ccall((:destroy_CMF_CS, testmod), Cvoid, (Ptr{Cvoid},), csf.cmf)
    end
end

function Base.finalize(csf::Field_C)
    destroy!(csf)
end

function load_config!(path::String)
    ccall((:load_config_julia, testmod), Cvoid, (Cstring,), path)
end

function (self::Field_C)(k::Vector{Float64}, w::Float64)::ComplexF64
    real_result::Ref{Float64} = Ref(Float64(0.0))
    imag_result::Ref{Float64} = Ref(Float64(0.0))
    ccall((:cmf_cs_call, testmod), Cvoid, (Ptr{Cvoid}, Ptr{Float64}, Float64, Cint, Ptr{Float64}, Ptr{Float64}), self.cmf, k, w, length(k), real_result, imag_result)
    return ComplexF64(real_result[], imag_result[])
end

function (self::Field_C)(w::Float64)::ComplexF64
    real_result::Ref{Float64} = Ref(Float64(0.0))
    imag_result::Ref{Float64} = Ref(Float64(0.0))
    ccall((:cmf_cs_call2, testmod), Cvoid, (Ptr{Cvoid}, Float64, Ptr{Float64}, Ptr{Float64}), self.cmf, w, real_result, imag_result)
    return ComplexF64(real_result[], imag_result[])
end

function save_field_to_file!(path::String, field)
    ccall((:cmf_save, testmod), Cvoid, (Cstring, Ptr{Cvoid}), path, field.cmf)
end

mutable struct Field_R
    cmf::Ptr{Cvoid}
    function Field_R(cmf::Ptr{Cvoid})
        new(cmf)
    end

    function Field_R(filename::String)
        cmf = ccall((:load_cmf_rs, testmod), Ptr{Cvoid}, (Cstring,), filename)
        if cmf == C_NULL
            error("Failed to load CMF_R from file: ", filename)
        end
        new(cmf)
    end
end

function (self::Field_R)(k::Vector{Float64}, w::Float64)::Float64
    real_result::Ref{Float64} = Ref(Float64(0.0))
    imag_result::Ref{Float64} = Ref(Float64(0.0))
    ccall((:cmf_rs_call, testmod), Cvoid, (Ptr{Cvoid}, Ptr{Float64}, Float64, Cint, Ptr{Float64}, Ptr{Float64}), self.cmf, k, w, length(k), real_result, imag_result)
    return Float64(real_result[])
end

function (self::Field_R)(w::Float64)::Float64
    real_result::Ref{Float64} = Ref(Float64(0.0))
    imag_result::Ref{Float64} = Ref(Float64(0.0))
    ccall((:cmf_rs_call2, testmod), Cvoid, (Ptr{Cvoid}, Float64, Ptr{Float64}, Ptr{Float64}), self.cmf, w, real_result, imag_result)
    return Float64(real_result[])
end

function Base.finalize(csf::Field_R)
    destroy!(csf)
end

function save_data(path::String, points, data, dimension, with_w, is_complex, is_vector)
    ccall((:save_to_file, testmod), Cvoid, (Cstring, Ptr{Float64}, Ptr{Float64}, Cint, Cint, Cint, Cint), path, points, data, dimension, with_w, is_complex, is_vector)
end

end # module Field

#using .Field

file = "/home/g/Research/Materials/Test/test_ckio_ir.dat"
#cmf2 = load_cmf_cs(file)
cmf2 = Field.Field_C(file)
if cmf2 == C_NULL
    error("Failed to load CMF_CS from file: ", file)
end
println(cmf2([-3.14, -3.14, -3.14], 1.0))
#println(cmf_cs_call(cmf2, [-3.14, -3.14, -3.14], 0.0))
#cmf3 = load_cmf_cs("/home/g/Research/Materials/Test/DOS.dat")
cmf3 = Field.Field_R("/home/g/Research/Materials/Test/test_DOS.dat")
if cmf3 == C_NULL
    error("Failed to load CMF_CS from file: DOS.dat")
end
#println(cmf_cs_call2(cmf3, 0.0))
#println(cmf3(0.0))
println(cmf3(1.0))

##path = "/home/g/Research/fcode/build/bin/input.cfg"
##load_config!(path)
##println(epsilon(1, [0.1, 0.2, 0.3]))
