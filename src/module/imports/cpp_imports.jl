#  g++ -shared -fPIC ../../build/CMakeFiles/ffirefly.x.dir/ffirefly/config/load/c_config.c.o ../../build/CMakeFiles/ffirefly.x.dir/ffirefly/objects/vec.cpp.o ../../build/CMakeFiles/ffirefly.x.dir/ffirefly/hamiltonian/band_structure.cpp.o ../../build/CMakeFiles/ffirefly.x.dir/ffirefly/config/load/cpp_config.cpp.o jmodtest.o -o jmodtest.so
module Imports

const libfly = abspath(@__FILE__)[1:end-51] * "build/lib/libfly.so"
export load_config!

# Begin Functions
export Bands, Field_C, Field_R, Vertex, epsilon

function epsilon(arg0::Int, arg1::Vector{Float64})
    newarg1 = Float32.(arg1)
    return ccall((:epsilon_export0, libfly), Float32, (Cint, Ptr{Float32}, Cint), arg0, newarg1, length(arg1))
end

mutable struct Bands
    ptr::Ptr{Cvoid}
end

function Bands()
    ptr = ccall((:Bands_export0, libfly), Ptr{Cvoid}, ())
    return Bands(ptr)
end

function (self::Bands)(arg0::Int, arg1::Vector{Float64})::Float32
    newarg1 = Float32.(arg1)
    lenarg1 = length(arg1)
    return ccall((:Bands_operator_export0, libfly), Float32, (Ptr{Cvoid}, Cint, Ptr{Float32}, Cint), self.ptr, arg0, newarg1, lenarg1)
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

function (self::Vertex)(arg0::Vector{Float64}, arg1=0.0, arg2::String="'up'", arg3::String="'up'")::ComplexF32
    newarg0 = Float32.(arg0)
    lenarg0 = length(arg0)
    newarg1 = Float32(arg1)
    real = Ref{Float32}()
    imag = Ref{Float32}()
    ccall((:Vertex_operator_export0, libfly), Nothing, (Ptr{Cvoid}, Ptr{Float32}, Cint, Float32, Cstring, Cstring, Ptr{Float32}, Ptr{Float32}), self.ptr, newarg0, lenarg0, newarg1, arg2, arg3, real, imag)
    return complex(real[], imag[])
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

function (self::Field_R)(arg0)::Float32
    newarg0 = Float32(arg0)
    return ccall((:Field_R_operator_export0, libfly), Float32, (Ptr{Cvoid}, Float32), self.ptr, newarg0)
end

function (self::Field_R)(arg0::Int, arg1)::Float32
    newarg1 = Float32(arg1)
    return ccall((:Field_R_operator_export1, libfly), Float32, (Ptr{Cvoid}, Cint, Float32), self.ptr, arg0, newarg1)
end

function (self::Field_R)(arg0::Vector{Float64}, arg1=0.0)::Float32
    newarg0 = Float32.(arg0)
    lenarg0 = length(arg0)
    newarg1 = Float32(arg1)
    return ccall((:Field_R_operator_export2, libfly), Float32, (Ptr{Cvoid}, Ptr{Float32}, Cint, Float32), self.ptr, newarg0, lenarg0, newarg1)
end

function (self::Field_R)(arg0::Int, arg1::Vector{Float64}, arg2=0.0)::Float32
    newarg1 = Float32.(arg1)
    lenarg1 = length(arg1)
    newarg2 = Float32(arg2)
    return ccall((:Field_R_operator_export3, libfly), Float32, (Ptr{Cvoid}, Cint, Ptr{Float32}, Cint, Float32), self.ptr, arg0, newarg1, lenarg1, newarg2)
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

function (self::Field_C)(arg0)::ComplexF32
    newarg0 = Float32(arg0)
    real = Ref{Float32}()
    imag = Ref{Float32}()
    ccall((:Field_C_operator_export0, libfly), Nothing, (Ptr{Cvoid}, Float32, Ptr{Float32}, Ptr{Float32}), self.ptr, newarg0, real, imag)
    return complex(real[], imag[])
end

function (self::Field_C)(arg0::Int, arg1)::ComplexF32
    newarg1 = Float32(arg1)
    real = Ref{Float32}()
    imag = Ref{Float32}()
    ccall((:Field_C_operator_export1, libfly), Nothing, (Ptr{Cvoid}, Cint, Float32, Ptr{Float32}, Ptr{Float32}), self.ptr, arg0, newarg1, real, imag)
    return complex(real[], imag[])
end

function (self::Field_C)(arg0::Vector{Float64}, arg1=0.0)::ComplexF32
    newarg0 = Float32.(arg0)
    lenarg0 = length(arg0)
    newarg1 = Float32(arg1)
    real = Ref{Float32}()
    imag = Ref{Float32}()
    ccall((:Field_C_operator_export2, libfly), Nothing, (Ptr{Cvoid}, Ptr{Float32}, Cint, Float32, Ptr{Float32}, Ptr{Float32}), self.ptr, newarg0, lenarg0, newarg1, real, imag)
    return complex(real[], imag[])
end

function (self::Field_C)(arg0::Int, arg1::Vector{Float64}, arg2=0.0)::ComplexF32
    newarg1 = Float32.(arg1)
    lenarg1 = length(arg1)
    newarg2 = Float32(arg2)
    real = Ref{Float32}()
    imag = Ref{Float32}()
    ccall((:Field_C_operator_export3, libfly), Nothing, (Ptr{Cvoid}, Cint, Ptr{Float32}, Cint, Float32, Ptr{Float32}, Ptr{Float32}), self.ptr, arg0, newarg1, lenarg1, newarg2, real, imag)
    return complex(real[], imag[])
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
