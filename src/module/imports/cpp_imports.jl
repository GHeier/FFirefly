#  g++ -shared -fPIC ../../build/CMakeFiles/ffirefly.x.dir/ffirefly/config/load/c_config.c.o ../../build/CMakeFiles/ffirefly.x.dir/ffirefly/objects/vec.cpp.o ../../build/CMakeFiles/ffirefly.x.dir/ffirefly/hamiltonian/band_structure.cpp.o ../../build/CMakeFiles/ffirefly.x.dir/ffirefly/config/load/cpp_config.cpp.o jmodtest.o -o jmodtest.so


module Imports

export epsilon,
       load_config!,
       Field_C,
       Field_R,
       destroy!,
       save_field_to_file!,
       save_data

const libfly = abspath(@__FILE__)[1:end-51] * "build/lib/libfly.so"

function epsilon(n::Int, k::Vector{Float64})
    newk = Float32.(k)
    return ccall((:epsilon_export0, libfly), Float32, (Int, Ptr{Float64}, Int), n, newk, length(k))
end

function load_config!(path::String)
    ccall((:load_config_export0, libfly), Cvoid, (Cstring,), path)
end

mutable struct Field_C
    ptr::Ptr{Cvoid}
end

function Field_C()
    ptr = ccall((:Field_C_export0, libfly), Ptr{Cvoid}, ())
  return Field_C(ptr)
end

function Field_C(filename::String)
    ptr = ccall((:Field_C_export2, libfly), Ptr{Cvoid}, (Cstring,), filename)
    return Field_C(ptr)
end

function (self::Field_C)(w)::ComplexF32
    real_result::Ref{Float32} = Ref(Float32(0.0))
    imag_result::Ref{Float32} = Ref(Float32(0.0))
    neww = Float32(w)
    ccall((:Field_C_operator_export0, libfly), Cvoid, (Ptr{Cvoid}, Cfloat, Ptr{Cfloat}, Ptr{Cfloat},), self.ptr, neww, real_result, imag_result)
    return ComplexF32(real_result[], imag_result[])
  end

function (self::Field_C)(n, w)::ComplexF32
    real_result::Ref{Float32} = Ref(Float32(0.0))
    imag_result::Ref{Float32} = Ref(Float32(0.0))
    neww = Float32(w)
    ccall((:Field_C_operator_export1, libfly), Cvoid, (Ptr{Cvoid}, Cint, Cfloat, Ptr{Cfloat}, Ptr{Cfloat},), self.ptr, n, neww, real_result, imag_result)
    return ComplexF32(real_result[], imag_result[])
  end

function (self::Field_C)(k::Vector{Float64}, w=0f0)::ComplexF32
    real_result::Ref{Float32} = Ref(Float32(0.0))
    imag_result::Ref{Float32} = Ref(Float32(0.0))
    newk::Vector{Float32} = Float32.(k)
    len = length(newk)
    neww = Float32(w)
    ccall((:Field_C_operator_export2, libfly), Cvoid, (Ptr{Cvoid}, Ptr{Float32}, Cint, Cfloat,Ptr{Cfloat}, Ptr{Cfloat},), self.ptr, newk, len, neww, real_result, imag_result)
    return ComplexF32(real_result[], imag_result[])
  end

function (self::Field_C)(n, k::Vector{Float64}, w=0f0)::ComplexF32
    real_result::Ref{Float32} = Ref(Float32(0.0))
    imag_result::Ref{Float32} = Ref(Float32(0.0))
    newk::Vector{Float32} = Float32.(k)
    len = length(newk)
    neww = Float32(w)
    ccall((:Field_C_operator_export3, libfly), Cvoid, (Ptr{Cvoid}, Ptr{Float32}, Cint, Cfloat,Ptr{Cfloat}, Ptr{Cfloat},), self.ptr, newk, len, neww, real_result, imag_result)
    return ComplexF32(real_result[], imag_result[])
  end

function destroy!(csf::Field_C)
    ccall((:destroy_Field_C, libfly), Cvoid, (Ptr{Cvoid},), csf.cmf)
end

function Base.finalize(csf::Field_C)
    destroy!(csf)
end

function save_field_to_file!(path::String, field)
    ccall((:cmf_save, libfly), Cvoid, (Cstring, Ptr{Cvoid}), path, field.cmf)
end

mutable struct Field_R
    ptr::Ptr{Cvoid}
end

function Field_R()
  ptr = ccall((:Field_R_export0, libfly), Ptr{Cvoid}, ())
  return Field_R(ptr)
end

function Field_R(filename::String)
    ptr = ccall((:Field_R_export2, libfly), Ptr{Cvoid}, (Cstring,), filename)
    return Field_R(ptr)
end

function (self::Field_R)(w)::Float32
    neww = Float32(w)
    return ccall((:Field_R_operator_export0, libfly), Cfloat, (Ptr{Cvoid}, Float32, ), self.ptr, neww)
end

function (self::Field_R)(n, w)::Float32
    neww = Float32(w)
    return ccall((:Field_R_operator_export1, libfly), Cfloat, (Ptr{Cvoid}, Cint, Float32, ), self.ptr, n, neww)
end

function (self::Field_R)(k::Vector{Float64}, w =0f0)::Float32
    neww = Float32(w)
    newk::Vector{Float32} = Float32.(k)
    len = length(newk)
    return ccall((:Field_R_operator_export2, libfly), Cfloat, (Ptr{Cvoid}, Ptr{Cfloat}, Cint, Float32, ), self.ptr, newk, len, neww)
end

function (self::Field_R)(n, k::Vector{Float64}, w=0f0)::Float32
    neww = Float32(w)
    newk::Vector{Float32} = Float32.(k)
    len = length(k)
    return ccall((:Field_R_operator_export3, libfly), Cfloat, (Ptr{Cvoid}, Cint, Ptr{Cfloat}, Cint, Float32, ), self.ptr, n, newk, len, neww)
end

function Base.finalize(csf::Field_R)
    destroy!(csf)
end

function save_data(path::String, points, data, dimension, with_w, is_complex, is_vector)
    ccall((:save_data, libfly), Cvoid, (Cstring, Ptr{Float64}, Ptr{Float64}, Cint, Cint, Cint, Cint), path, points, data, dimension, with_w, is_complex, is_vector)
end

end # module Field

#using .Imports

#testcmf = Ffirefly.Field_R("../../../sample_bands.dat")
#println(testcmf(1, [-1.0, -0.6]))
#
#path = "/home/g/Research/FFirefly/build/bin/input.cfg"
#load_config!(path)
#println(epsilon(1, [0.1, 0.2, 0.3]))
