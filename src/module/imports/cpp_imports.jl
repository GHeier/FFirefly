#  g++ -shared -fPIC ../../build/CMakeFiles/ffirefly.x.dir/ffirefly/config/load/c_config.c.o ../../build/CMakeFiles/ffirefly.x.dir/ffirefly/objects/vec.cpp.o ../../build/CMakeFiles/ffirefly.x.dir/ffirefly/hamiltonian/band_structure.cpp.o ../../build/CMakeFiles/ffirefly.x.dir/ffirefly/config/load/cpp_config.cpp.o jmodtest.o -o jmodtest.so
module Imports

const libfly = abspath(@__FILE__)[1:end-51] * "build/lib/libfly.so"
export load_config!, Vec, Surface, get_faces
export epsilon,
       vk,
       norm,
       Vertex,
       Self_Energy,
       Field_C,
       Field_R,
       destroy!,
       save_field_to_file!,
       save_data,
       data_save!

# Struct for Vec
struct RawVec
    x::Float32
    y::Float32
    z::Float32
    w::Float32
    area::Float32
    dimension::Int32
    n::Int32
end


# Vec class
mutable struct Vec
    ptr::Ptr{Cvoid}
    x::Float32
    y::Float32
    z::Float32
    w::Float32
    area::Float32
    dimension::Int32
    n::Int32

    function Vec(args...)
        ptr = C_NULL
        if length(args) == 0
            ptr = ccall((:Vec_export0, libfly), Ptr{Cvoid}, ())
        elseif length(args) >= 1 && (args[1] isa Float32 || args[1] isa Float64)
            x = Float64(args[1])
            y = length(args) > 1 ? Float64(args[2]) : 0
            z = length(args) > 2 ? Float64(args[3]) : 0
            w = length(args) > 3 ? Float64(args[4]) : 0
            area = length(args) > 4 ? Float64(args[5]) : 0
            dimension = length(args) > 5 ? args[6] : 3
            n = length(args) > 6 ? args[7] : 0
            ptr = ccall((:Vec_export1, libfly), Ptr{Cvoid},
                        (Cfloat, Cfloat, Cfloat, Cfloat, Cfloat, Cint, Cint),
                        x, y, z, w, area, dimension, n)
        end

        if ptr == C_NULL
            error("Failed to initialize Vec")
        end

        x = ccall((:Vec_x_export0, libfly), Cfloat, (Ptr{Cvoid},), ptr)
        y = ccall((:Vec_y_export0, libfly), Cfloat, (Ptr{Cvoid},), ptr)
        z = ccall((:Vec_z_export0, libfly), Cfloat, (Ptr{Cvoid},), ptr)
        w = ccall((:Vec_w_export0, libfly), Cfloat, (Ptr{Cvoid},), ptr)
        area = ccall((:Vec_area_export0, libfly), Cfloat, (Ptr{Cvoid},), ptr)
        dimension = ccall((:Vec_dimension_export0, libfly), Cint, (Ptr{Cvoid},), ptr)
        n = ccall((:Vec_n_export0, libfly), Cint, (Ptr{Cvoid},), ptr)

        new(ptr, x, y, z, w, area, dimension, n)
    end

    function finalize(v::Vec)
        try
            ccall((:destroy, libfly), Cvoid, (Ptr{Cvoid},), v.ptr)
        catch
        end
    end
end

Base.:+(a::Vec, b::Vec) = Vec(a.x + b.x, a.y + b.y, a.z + b.z)
Base.:-(a::Vec, b::Vec) = Vec(a.x - b.x, a.y - b.y, a.z - b.z)
Base.:*(a, b::Vec) = Vec(b.x * a, b.y * a, b.z * a)
Base.:/(a, b::Vec) = Vec(b.x / a, b.y / a, b.z / a)
Base.:*(b::Vec, a) = Vec(b.x * a, b.y * a, b.z * a)
Base.:/(b::Vec, a) = Vec(b.x / a, b.y / a, b.z / a)

# Utility functions
function string_to_vec(s::String)
    return ccall((:string_to_vec_export0, libfly), Ptr{Cvoid}, (Cstring,), s)
end

function unpack_string(s::String)
    return ccall((:unpack_string_export0, libfly), Ptr{Cvoid}, (Cstring,), s)
end

function vec_to_string(vec::Vector{Float32})
    return unsafe_string(ccall((:vec_to_string_export0, libfly), Cstring,
                               (Ptr{Cfloat}, Cint), vec, length(vec)))
end

function round_(n::Int)
    return ccall((:round_export0, libfly), Ptr{Cvoid}, (Cint,), n)
end

function norm(v::Vec)
    val = RawVec(v.x, v.y, v.z, v.w, v.area, v.dimension, v.n)
    return ccall((:norm_export0, libfly), Cfloat, (Ref{RawVec},), val)
end

mutable struct Surface
    handle::Ptr{Cvoid}
    faces::Vector{Vec}
end

const _userfunc_registry = IdDict{Ptr{Cvoid}, Function}()
const current_callback_key = Ref{Ptr{Cvoid}}(C_NULL)

function _dispatch_callback(ptr::Ptr{RawVec})::Float32
    k = unsafe_load(ptr)
    func = _userfunc_registry[current_callback_key[]]
    return func(k)
end

const _trampoline = @cfunction(_dispatch_callback, Float32, (Ptr{RawVec},))
function Surface(userfunc::Function, s_val)
    s_val = Float32(s_val)
    key = Base.unsafe_convert(Ptr{Cvoid}, Ref(userfunc))  # unique key
    _userfunc_registry[key] = userfunc
    current_callback_key[] = key

    handle = ccall((:Surface_export0, libfly), Ptr{Cvoid},
                (Ptr{Cvoid}, Cfloat), _trampoline, s_val)
    faces = get_faces1(handle)

    return Surface(handle, faces)
end

function get_faces(surf::Surface)::Vector{Vector{Float32}}
    # Step 1: Get number of faces
    n_faces = ccall((:Surface_num_faces_export0, libfly), Cint,
                    (Ptr{Cvoid},), surf.handle)

    if n_faces <= 0
        return []
    end

    # Step 2: Prepare buffers
    lens = Vector{Cint}(undef, n_faces)
    total_len = 3 * n_faces  # Adjust if needed based on your data shape
    buf = Vector{Cfloat}(undef, total_len)
    n_faces_ref = Ref{Cint}(n_faces)

    # Step 3: Call C++ function
    ccall((:Surface_var_faces_export0, libfly), Cvoid,
          (Ptr{Cvoid}, Ptr{Cfloat}, Ptr{Cint}, Ptr{Cint}),
          surf.handle, buf, lens, n_faces_ref)
    # Step 4: Reconstruct nested vector
    result = Vector{Vector{Float32}}()
    offset = 0
    for i in 1:n_faces
        len = lens[i]
        push!(result, buf[offset+1 : offset+len])
        offset += len
    end

    return result
end

function get_faces1(handle::Ptr{Cvoid})::Vector{Vec}
    # Step 1: Get number of faces
    n_faces = ccall((:Surface_num_faces_export0, libfly), Cint,
                    (Ptr{Cvoid},), handle)

    if n_faces <= 0
        return []
    end

    # Step 2: Prepare buffers
    buf = Vector{RawVec}(undef, n_faces)
    # Step 3: Call C++ function
    ccall((:Surface_var_faces_export1, libfly), Cvoid,
          (Ptr{Cvoid}, Ptr{RawVec}),
          handle, buf)
    result = Vector{Vec}(undef, n_faces)
    for i in 1:n_faces
        b = buf[i]
        result[i] = Vec(b.x, b.y, b.z, b.w, b.area, b.dimension, b.n)
    end
    return result
end


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

function (self::Bands)(arg0::Int, arg1::RawVec)::Float32
    lenarg1 = arg1.dimension
    newarg1 = zeros(Float32, lenarg1)
    newarg1[1] = arg1.x
    if lenarg1 > 1 
        newarg1[2] = arg1.y
        if lenarg1 > 2 
            newarg1[3] = arg1.z
        end 
    end
    return ccall((:Bands_operator_export0, libfly), Float32, (Ptr{Cvoid}, Cint, Ptr{Float32}, Cint), self.ptr, arg0, newarg1, lenarg1)
end

function vk(n, k::Vec, band::Bands)::Vec
    nint = Int32(n)
    newk = RawVec(k.x, k.y, k.z, k.w, k.area, k.dimension, k.n)
    result = ccall((:vk_export0, libfly), Ptr{RawVec}, (Cint, Ref{RawVec}, Ptr{Cvoid}), nint, newk, band.ptr) 
    q = unsafe_load(result)
    return Vec(q.x, q.y, q.z, q.w, q.area, q.dimension, q.n)
end

function Base.finalize(obj::Bands)
    destroy!(obj)
end

mutable struct Self_Energy
    ptr::Ptr{Cvoid}
end

function Self_Energy()
    ptr = ccall((:Self_Energy_export0, libfly), Ptr{Cvoid}, ())
    return Self_Energy(ptr)
end

function (self::Self_Energy)(arg0::Vector{Float64}, arg1::Float64)::ComplexF32
    real_result = Ref{Cfloat}(0.0f0)
    imag_result = Ref{Cfloat}(0.0f0)
    newarg0 = Float32.(arg0)
    lenarg0 = length(newarg0)  # use the converted arg
    ccall((:Self_Energy_operator_export0, libfly), Cvoid,
          (Ptr{Cvoid}, Ptr{Cfloat}, Cint, Cfloat, Ptr{Cfloat}, Ptr{Cfloat}),
          self.ptr, newarg0, lenarg0, Float32(arg1), real_result, imag_result)
    return ComplexF32(real_result[], imag_result[])
end

function (self::Self_Energy)(arg0::Vec, arg1::Float64)::ComplexF32
    real_result = Ref{Cfloat}(0.0f0)
    imag_result = Ref{Cfloat}(0.0f0)
    newarg0::Vector{Float32} = [arg0.x, arg0.y, arg0.z]
    lenarg0 = length(newarg0)  # use the converted arg
    ccall((:Self_Energy_operator_export0, libfly), Cvoid,
          (Ptr{Cvoid}, Ptr{Cfloat}, Cint, Cfloat, Ptr{Cfloat}, Ptr{Cfloat}),
          self.ptr, newarg0, lenarg0, Float32(arg1), real_result, imag_result)
    return ComplexF32(real_result[], imag_result[])
end

function Base.finalize(obj::Self_Energy)
    destroy!(obj)
end

mutable struct Vertex
    ptr::Ptr{Cvoid}
end

function Vertex()
    ptr = ccall((:Vertex_export0, libfly), Ptr{Cvoid}, ())
    return Vertex(ptr)
end

function (self::Vertex)(k::Vector{Float64}, w=0f0)::ComplexF32
    real_result::Ref{Float32} = Ref(Float32(0.0))
    imag_result::Ref{Float32} = Ref(Float32(0.0))
    newk::Vector{Float32} = Float32.(k)
    len = length(newk)
    neww = Float32(w)
    ccall((:Vertex_operator_export0, libfly), Cvoid, (Ptr{Cvoid}, Ptr{Float32}, Cint, Cfloat,Ptr{Cfloat}, Ptr{Cfloat},), self.ptr, newk, len, neww, real_result, imag_result)
    return ComplexF32(real_result[], imag_result[])
end

function (self::Vertex)(k::Vec, w=0f0)::ComplexF32
    real_result::Ref{Float32} = Ref(Float32(0.0))
    imag_result::Ref{Float32} = Ref(Float32(0.0))
    newk::Vector{Float32} = [k.x, k.y, k.z]
    len = length(newk)
    neww = Float32(w)
    ccall((:Vertex_operator_export0, libfly), Cvoid, (Ptr{Cvoid}, Ptr{Float32}, Cint, Cfloat,Ptr{Cfloat}, Ptr{Cfloat},), self.ptr, newk, len, neww, real_result, imag_result)
    return ComplexF32(real_result[], imag_result[])
end

function Base.finalize(obj::Vertex)
    destroy!(obj)
end

mutable struct Self_Energy
    ptr::Ptr{Cvoid}
end

function Self_Energy()
    ptr = ccall((:Self_Energy_export0, libfly), Ptr{Cvoid}, ())
    return Self_Energy(ptr)
end

function (self::Self_Energy)(k::Vector{Float64}, w=0f0)::ComplexF32
    real_result::Ref{Float32} = Ref(Float32(0.0))
    imag_result::Ref{Float32} = Ref(Float32(0.0))
    newk::Vector{Float32} = Float32.(k)
    len = length(newk)
    neww = Float32(w)
    ccall((:Self_Energy_operator_export0, libfly), Cvoid, (Ptr{Cvoid}, Ptr{Float32}, Cint, Cfloat,Ptr{Cfloat}, Ptr{Cfloat},), self.ptr, newk, len, neww, real_result, imag_result)
    return ComplexF32(real_result[], imag_result[])
end

function (self::Self_Energy)(k::Vec, w=0f0)::ComplexF32
    real_result::Ref{Float32} = Ref(Float32(0.0))
    imag_result::Ref{Float32} = Ref(Float32(0.0))
    newk::Vector{Float32} = [k.x, k.y, k.z]
    len = length(newk)
    neww = Float32(w)
    ccall((:Self_Energy_operator_export0, libfly), Cvoid, (Ptr{Cvoid}, Ptr{Float32}, Cint, Cfloat,Ptr{Cfloat}, Ptr{Cfloat},), self.ptr, newk, len, neww, real_result, imag_result)
    return ComplexF32(real_result[], imag_result[])
end

function Base.finalize(obj::Self_Energy)
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
    ptr = ccall((:Field_R_export2, libfly), Ptr{Cvoid}, (Cstring,), filename)
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

function destroy!(csf::Field_C)
    ccall((:destroy_Field_C, libfly), Cvoid, (Ptr{Cvoid},), csf.cmf)
end

function Base.finalize(obj::Field_C)
    destroy!(obj)
end

# End Functions

function save_field_to_file!(path::String, field)
    ccall((:cmf_save, libfly), Cvoid, (Cstring, Ptr{Cvoid}), path, field.cmf)
end

function save_data(path::String, points, data, dimension, with_w, is_complex, is_vector)
    ccall((:save_data, libfly), Cvoid, (Cstring, Ptr{Float64}, Ptr{Float64}, Cint, Cint, Cint, Cint), path, points, data, dimension, with_w, is_complex, is_vector)
end

function data_save!(path::String, points, data, dimension, with_w, is_complex, is_vector)
    println("type of points: ", typeof(points))
    println("type of data: ", typeof(data))
    println("type of dimension: ", typeof(dimension))
    numpts = length(points)
    jpoints = points
    jdata = data
    ccall((:data_save_export0, libfly), Cvoid, (Cstring, Ptr{Float64}, Ptr{Float64}, Cint, Cint, Cint, Cint, Cint), path, jpoints, jdata, numpts, dimension, with_w, is_complex, is_vector)
end

function load_config!(path::String)
    ccall((:load_config_export0, libfly), Cvoid, (Cstring,), path)
end

end # module Imports

#using .Imports

#testcmf = Ffirefly.Field_R("../../../sample_bands.dat")
#println(testcmf(1, [-1.0, -0.6]))
#
#path = "/home/g/Research/FFirefly/build/bin/input.cfg"
#load_config!(path)
#println(epsilon(1, [0.1, 0.2, 0.3]))
