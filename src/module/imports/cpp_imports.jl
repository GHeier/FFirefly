#  g++ -shared -fPIC ../../build/CMakeFiles/ffirefly.x.dir/ffirefly/config/load/c_config.c.o ../../build/CMakeFiles/ffirefly.x.dir/ffirefly/objects/vec.cpp.o ../../build/CMakeFiles/ffirefly.x.dir/ffirefly/hamiltonian/band_structure.cpp.o ../../build/CMakeFiles/ffirefly.x.dir/ffirefly/config/load/cpp_config.cpp.o jmodtest.o -o jmodtest.so
module Imports

const libfly = abspath(@__FILE__)[1:end-51] * "build/lib/libfly.so"
export load_config!

# Begin Functions
export Bands, Field_C, Field_C, Field_R, Field_R, Vertex, epsilon

function epsilon(arg0::Int, arg1::Vector{Float64})::Cfloat
    newarg1 = Float32.(arg1)
    return ccall((:epsilon_export0, libfly), Cfloat, (Cint, Ptr{Cfloat}, Cint), arg0, newarg1, length(arg1))
end

mutable struct Bands
    ptr::Ptr{Cvoid}
end

function Bands()
    ptr = ccall((:Bands_export0, libfly), Ptr{Cvoid}, ())
    return Bands(ptr)
end

function Base.finalize(obj::Bands)
    destroy!(obj)
end

mutable struct Vertex
    ptr::Ptr{Cvoid}
end

function Vertex()
    ptr = ccall((:Vertex_export0, libfly), Ptr{Cvoid}, ())
    return Vertex(ptr)
end

function Base.finalize(obj::Vertex)
    destroy!(obj)
end

mutable struct Field_R
    ptr::Ptr{Cvoid}
end

function Field_R()
    ptr = ccall((:Field_R_export0, libfly), Ptr{Cvoid}, ())
    return Field_R(ptr)
end

function Field_R(filename::String)
    ptr = ccall((:Field_R_export1, libfly), Ptr{Cvoid}, (Cstring,), filename)
    return Field_R(ptr)
end

function (self::Field_R)(arg0::Float32)::Cfloat
    newarg0 = Cfloat(arg0)
    return ccall((:Field_R_operator_export0, libfly), Cfloat, (Ptr{Cvoid}, Cfloat), self.ptr, newarg0)
end

function (self::Field_R)(arg0::Int, arg1::Float32)::Cfloat
    newarg1 = Cfloat(arg1)
    return ccall((:Field_R_operator_export1, libfly), Cfloat, (Ptr{Cvoid}, Cint, Cfloat), self.ptr, arg0, newarg1)
end

function (self::Field_R)(arg0::Vector{Float64}, arg1::Float32 = 0.0)::Cfloat
    newarg0 = Float32.(arg0)
    len = length(arg0)
    newarg1 = Cfloat(arg1)
    return ccall((:Field_R_operator_export2, libfly), Cfloat, (Ptr{Cvoid}, Ptr{Cfloat}, Cint, Cfloat), self.ptr, newarg0, len, newarg1)
end

function (self::Field_R)(arg0::Int, arg1::Vector{Float64}, arg2::Float32 = 0.0)::Cfloat
    newarg1 = Float32.(arg1)
    len = length(arg1)
    newarg2 = Cfloat(arg2)
    return ccall((:Field_R_operator_export3, libfly), Cfloat, (Ptr{Cvoid}, Cint, Ptr{Cfloat}, Cint, Cfloat), self.ptr, arg0, newarg1, len, newarg2)
end

function Base.finalize(obj::Field_R)
    destroy!(obj)
end

mutable struct Field_C
    ptr::Ptr{Cvoid}
end

function Field_C()
    ptr = ccall((:Field_C_export0, libfly), Ptr{Cvoid}, ())
    return Field_C(ptr)
end

function Field_C(filename::String)
    ptr = ccall((:Field_C_export1, libfly), Ptr{Cvoid}, (Cstring,), filename)
    return Field_C(ptr)
end

function (self::Field_C)(arg0::Float32)::ComplexF32
    newarg0 = Cfloat(arg0)
    re = Ref{Cfloat}()
    im = Ref{Cfloat}()
    ccall((:Field_C_operator_export0, libfly), Nothing, (Ref{Cfloat}, Ref{Cfloat}, Ptr{Cvoid}, Cfloat), re, im, self.ptr, newarg0)
    return ComplexF32(re[], im[])
end

function (self::Field_C)(arg0::Int, arg1::Float32)::ComplexF32
    newarg1 = Cfloat(arg1)
    re = Ref{Cfloat}()
    im = Ref{Cfloat}()
    ccall((:Field_C_operator_export1, libfly), Nothing, (Ref{Cfloat}, Ref{Cfloat}, Ptr{Cvoid}, Cint, Cfloat), re, im, self.ptr, arg0, newarg1)
    return ComplexF32(re[], im[])
end

function (self::Field_C)(arg0::Vector{Float64}, arg1::Float32 = 0.0)::ComplexF32
    newarg0 = Float32.(arg0)
    len = length(arg0)
    newarg1 = Cfloat(arg1)
    re = Ref{Cfloat}()
    im = Ref{Cfloat}()
    ccall((:Field_C_operator_export2, libfly), Nothing, (Ref{Cfloat}, Ref{Cfloat}, Ptr{Cvoid}, Ptr{Cfloat}, Cint, Cfloat), re, im, self.ptr, newarg0, len, newarg1)
    return ComplexF32(re[], im[])
end

function (self::Field_C)(arg0::Int, arg1::Vector{Float64}, arg2::Float32 = 0.0)::ComplexF32
    newarg1 = Float32.(arg1)
    len = length(arg1)
    newarg2 = Cfloat(arg2)
    re = Ref{Cfloat}()
    im = Ref{Cfloat}()
    ccall((:Field_C_operator_export3, libfly), Nothing, (Ref{Cfloat}, Ref{Cfloat}, Ptr{Cvoid}, Cint, Ptr{Cfloat}, Cint, Cfloat), re, im, self.ptr, arg0, newarg1, len, newarg2)
    return ComplexF32(re[], im[])
end

function Base.finalize(obj::Field_C)
    destroy!(obj)
end

# End Functions

function save_data(path::String, points, data, dimension, with_w, is_complex, is_vector)
    ccall((:save_data, libfly), Cvoid, (Cstring, Ptr{Float64}, Ptr{Float64}, Cint, Cint, Cint, Cint), path, points, data, dimension, with_w, is_complex, is_vector)
end

function load_config!(path::String)
    ccall((:load_config_export0, libfly), Cvoid, (Cstring,), path)
end

end # module Field

#using .Imports

#testcmf = Ffirefly.Field_R("../../../sample_bands.dat")
#println(testcmf(1, [-1.0, -0.6]))
#
#path = "/home/g/Research/FFirefly/build/bin/input.cfg"
#load_config!(path)
#println(epsilon(1, [0.1, 0.2, 0.3]))
