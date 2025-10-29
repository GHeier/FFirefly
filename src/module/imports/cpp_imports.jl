#  g++ -shared -fPIC ../../build/CMakeFiles/ffirefly.x.dir/ffirefly/config/load/c_config.c.o ../../build/CMakeFiles/ffirefly.x.dir/ffirefly/objects/vec.cpp.o ../../build/CMakeFiles/ffirefly.x.dir/ffirefly/hamiltonian/band_structure.cpp.o ../../build/CMakeFiles/ffirefly.x.dir/ffirefly/config/load/cpp_config.cpp.o jmodtest.o -o jmodtest.so
module Imports

# Get the project root directory (FFirefly/)
const libfly = split(abspath(@__FILE__), "FFirefly")[1] * "FFirefly/build/lib/libfly.so"
#const libfly = joinpath(dirname(dirname(dirname(dirname(abspath(@__FILE__))))), "build", "lib", "libfly.so")
export load_config!, Vec, Surface, get_faces
export epsilon,
       norm,
       Bands,
       Vertex,
       Self_Energy,
       Hamiltonian,
       Field_C,
       Field_R,
       Field_RM,
       Field_CM,
       destroy!,
       save_field_to_file!,
       save_data!,
       data_save!,
       save_field!,
       save_data_scalar,
       save_data_vector,
       save_data_matrix

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
export Field_C, Field_R, Field_RM, Field_CM, Vertex, epsilon

function epsilon(arg0::Int, arg1::Vector{Float64})
    newarg1 = Float32.(arg1)
    return ccall((:epsilon_export0, libfly), Float32, (Cint, Ptr{Float32}, Cint), arg0, newarg1, length(arg1))
end

mutable struct Bands
    ptr::Ptr{Cvoid}
end

function Bands()
    println("Julia Bands Constructor called")
    ptr = ccall((:Bands_export0, libfly), Ptr{Cvoid}, ())
    println("Julia Bands Constructor finished")
    return Bands(ptr)
end

function (self::Bands)(arg0::Int, arg1::Vector{Float64})::Float32
    newarg1 = Float32.(arg1)
    lenarg1 = length(arg1)
    return ccall((:Bands_operator_export0, libfly), Float32, (Ptr{Cvoid}, Cint, Ptr{Float32}, Cint), self.ptr, arg0, newarg1, lenarg1)
end

function (self::Bands)(arg0::Int, arg1::Vec)::Float32
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

function (self::Bands)(arg1::Vector{Float64})::Float32
    newarg1 = Float32.(arg1)
    lenarg1 = length(arg1)
    return ccall((:Bands_operator_export1, libfly), Float32, (Ptr{Cvoid}, Ptr{Float32}, Cint), self.ptr, newarg1, lenarg1)
end

function (self::Bands)(arg1::Vec)::Float32
    lenarg1 = arg1.dimension
    newarg1 = zeros(Float32, lenarg1)
    newarg1[1] = arg1.x
    if lenarg1 > 1
        newarg1[2] = arg1.y
        if lenarg1 > 2
            newarg1[3] = arg1.z
        end
    end
    return ccall((:Bands_operator_export1, libfly), Float32, (Ptr{Cvoid}, Ptr{Float32}, Cint), self.ptr, newarg1, lenarg1)
end

function Base.finalize(obj::Bands)
    ccall((:destroy_Bands, libfly), Cvoid, (Ptr{Cvoid},), obj.ptr)
end

mutable struct Self_Energy
    ptr::Ptr{Cvoid}
end

function Self_Energy()
    ptr = ccall((:Self_Energy_export0, libfly), Ptr{Cvoid}, ())
    return Self_Energy(ptr)
end

function (self::Self_Energy)(arg0::Vector{Float64}, arg1)::ComplexF32
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

mutable struct Field_R
    ptr::Ptr{Cvoid}
    dimension::Int
    mesh::Vector{Int}
    domain::Matrix{Float32}
    w_points::Vector{Float32}
end

function Field_R()
    ptr = ccall((:Field_R_export0, libfly), Ptr{Cvoid}, ())
    obj = Field_R(ptr, 0, Int[], zeros(Float32, 0, 0), Float32[])
    _load_field_metadata!(obj)
    return obj
end

function Field_R(filename::String)
    ptr = ccall((:Field_R_export2, libfly), Ptr{Cvoid}, (Cstring,), filename)
    obj = Field_R(ptr, 0, Int[], zeros(Float32, 0, 0), Float32[])
    _load_field_metadata!(obj)
    return obj
end

function _load_field_metadata!(obj::Field_R)
    # Get dimension
    obj.dimension = ccall((:Field_R_get_dimension, libfly), Cint, (Ptr{Cvoid},), obj.ptr)

    # Get mesh
    mesh_size = ccall((:Field_R_get_mesh_size, libfly), Cint, (Ptr{Cvoid},), obj.ptr)
    if mesh_size > 0
        obj.mesh = zeros(Int32, mesh_size)
        ccall((:Field_R_get_mesh, libfly), Cvoid, (Ptr{Cvoid}, Ptr{Cint}), obj.ptr, obj.mesh)
        obj.mesh = Int.(obj.mesh)
    end

    # Get domain
    domain_rows = ccall((:Field_R_get_domain_rows, libfly), Cint, (Ptr{Cvoid},), obj.ptr)
    domain_cols = ccall((:Field_R_get_domain_cols, libfly), Cint, (Ptr{Cvoid},), obj.ptr)
    if domain_rows > 0 && domain_cols > 0
        domain_flat = zeros(Float32, domain_rows * domain_cols)
        ccall((:Field_R_get_domain, libfly), Cvoid, (Ptr{Cvoid}, Ptr{Float32}), obj.ptr, domain_flat)
        obj.domain = reshape(domain_flat, domain_cols, domain_rows)'
    end

    # Get w_points
    w_points_size = ccall((:Field_R_get_w_points_size, libfly), Cint, (Ptr{Cvoid},), obj.ptr)
    if w_points_size > 0
        obj.w_points = zeros(Float32, w_points_size)
        ccall((:Field_R_get_w_points, libfly), Cvoid, (Ptr{Cvoid}, Ptr{Float32}), obj.ptr, obj.w_points)
    end
end

function (self::Field_R)(arg0)::Float32
    newarg0 = Float32(arg0)
    return ccall((:Field_R_operator_export0, libfly), Float32, (Ptr{Cvoid}, Float32), self.ptr, newarg0)
end

# Multiple w-points evaluation - explicit function name to avoid dispatch ambiguity
function eval_w_list(self::Field_R, w_points::AbstractVector{<:Real})::Vector{Float32}
    num_w = length(w_points)
    if num_w == 0
        return Float32[]
    end

    w_array = Float32.(w_points)
    output = zeros(Float32, num_w)

    ccall((:Field_R_operator_export_w_list, libfly), Cvoid,
          (Ptr{Cvoid}, Ptr{Float32}, Cint, Ptr{Float32}),
          self.ptr, w_array, num_w, output)

    return output
end

# COMMENTED OUT - no corresponding export
#function (self::Field_R)(arg0::Int, arg1)::Float32
#    newarg1 = Float32(arg1)
#    return ccall((:Field_R_operator_export1, libfly), Float32, (Ptr{Cvoid}, Cint, Float32), self.ptr, arg0, newarg1)
#end

function (self::Field_R)(arg0::Vector{Float64}, arg1=0.0)::Float32
    newarg0 = Float32.(arg0)
    lenarg0 = length(arg0)
    newarg1 = Float32(arg1)
    return ccall((:Field_R_operator_export2, libfly), Float32, (Ptr{Cvoid}, Ptr{Float32}, Cint, Float32), self.ptr, newarg0, lenarg0, newarg1)
end

function (self::Field_R)(points::Vector{Vector{Float64}}, w=0.0)::Vector{Float32}
    num_points = length(points)
    if num_points == 0
        return Float32[]
    end

    len = length(points[1])
    points_flat = Float32[]
    for p in points
        append!(points_flat, Float32.(p))
    end

    output = zeros(Float32, num_points)
    ccall((:Field_R_operator_export_list, libfly), Cvoid,
          (Ptr{Cvoid}, Ptr{Float32}, Cint, Cint, Float32, Ptr{Float32}),
          self.ptr, points_flat, num_points, len, Float32(w), output)
    return output
end

# COMMENTED OUT - no corresponding export
#function (self::Field_R)(arg0::Int, arg1::Vector{Float64}, arg2=0.0)::Float32
#    newarg1 = Float32.(arg1)
#    lenarg1 = length(arg1)
#    newarg2 = Float32(arg2)
#    return ccall((:Field_R_operator_export3, libfly), Float32, (Ptr{Cvoid}, Cint, Ptr{Float32}, Cint, Float32), self.ptr, arg0, newarg1, lenarg1, newarg2)
#end

function Base.finalize(obj::Field_R)
    destroy!(obj)
end

function get_data(obj::Field_R)
    """Get data array reshaped in (w,k) format."""
    ptr = ccall((:Field_R_get_data, libfly), Ptr{Cvoid}, (Ptr{Cvoid},), obj.ptr)
    bd = _basedata_from_ptr(ptr)
    data_array = get_data(bd)

    # Reshape from (nk*nw) to (nw, nk)
    if bd.n_indices == 2
        # Matrix: (dim, dim, nk*nw) -> (dim, dim, nk, nw) -> (dim, dim, nw, nk)
        data_reshaped = reshape(data_array, bd.dim_indices, bd.dim_indices, bd.nk, bd.nw)
        return permutedims(data_reshaped, (1, 2, 4, 3))  # swap k and w axes
    else
        # Scalar: (nk*nw,) -> (nk, nw) -> (nw, nk)
        data_reshaped = reshape(data_array, bd.nk, bd.nw)
        return permutedims(data_reshaped, (2, 1))  # transpose to (nw, nk)
    end
end

mutable struct Field_C
    ptr::Ptr{Cvoid}
    dimension::Int
    mesh::Vector{Int}
    domain::Matrix{Float32}
    w_points::Vector{Float32}
end

function Field_C()
    ptr = ccall((:Field_C_export0, libfly), Ptr{Cvoid}, ())
    obj = Field_C(ptr, 0, Int[], zeros(Float32, 0, 0), Float32[])
    _load_field_metadata!(obj)
    return obj
end

function Field_C(filename::String)
    ptr = ccall((:Field_C_export2, libfly), Ptr{Cvoid}, (Cstring,), filename)
    obj = Field_C(ptr, 0, Int[], zeros(Float32, 0, 0), Float32[])
    _load_field_metadata!(obj)
    return obj
end

function _load_field_metadata!(obj::Field_C)
    # Get dimension
    obj.dimension = ccall((:Field_C_get_dimension, libfly), Cint, (Ptr{Cvoid},), obj.ptr)

    # Get mesh
    mesh_size = ccall((:Field_C_get_mesh_size, libfly), Cint, (Ptr{Cvoid},), obj.ptr)
    if mesh_size > 0
        obj.mesh = zeros(Int32, mesh_size)
        ccall((:Field_C_get_mesh, libfly), Cvoid, (Ptr{Cvoid}, Ptr{Cint}), obj.ptr, obj.mesh)
        obj.mesh = Int.(obj.mesh)
    end

    # Get domain
    domain_rows = ccall((:Field_C_get_domain_rows, libfly), Cint, (Ptr{Cvoid},), obj.ptr)
    domain_cols = ccall((:Field_C_get_domain_cols, libfly), Cint, (Ptr{Cvoid},), obj.ptr)
    if domain_rows > 0 && domain_cols > 0
        domain_flat = zeros(Float32, domain_rows * domain_cols)
        ccall((:Field_C_get_domain, libfly), Cvoid, (Ptr{Cvoid}, Ptr{Float32}), obj.ptr, domain_flat)
        obj.domain = reshape(domain_flat, domain_cols, domain_rows)'
    end

    # Get w_points
    w_points_size = ccall((:Field_C_get_w_points_size, libfly), Cint, (Ptr{Cvoid},), obj.ptr)
    if w_points_size > 0
        obj.w_points = zeros(Float32, w_points_size)
        ccall((:Field_C_get_w_points, libfly), Cvoid, (Ptr{Cvoid}, Ptr{Float32}), obj.ptr, obj.w_points)
    end
end

function (self::Field_C)(arg0)::ComplexF32
    newarg0 = Float32(arg0)
    real = Ref{Float32}()
    imag = Ref{Float32}()
    ccall((:Field_C_operator_export0, libfly), Nothing, (Ptr{Cvoid}, Float32, Ptr{Float32}, Ptr{Float32}), self.ptr, newarg0, real, imag)
    return complex(real[], imag[])
end

# Multiple w-points evaluation - explicit function name to avoid dispatch ambiguity
function eval_w_list(self::Field_C, w_points::AbstractVector{<:Real})::Vector{ComplexF32}
    num_w = length(w_points)
    if num_w == 0
        return ComplexF32[]
    end

    w_array = Float32.(w_points)
    real_output = zeros(Float32, num_w)
    imag_output = zeros(Float32, num_w)

    ccall((:Field_C_operator_export_w_list, libfly), Cvoid,
          (Ptr{Cvoid}, Ptr{Float32}, Cint, Ptr{Float32}, Ptr{Float32}),
          self.ptr, w_array, num_w, real_output, imag_output)

    return complex.(real_output, imag_output)
end

export eval_w_list

# COMMENTED OUT - no corresponding export
#function (self::Field_C)(arg0::Int, arg1)::ComplexF32
#    newarg1 = Float32(arg1)
#    real = Ref{Float32}()
#    imag = Ref{Float32}()
#    ccall((:Field_C_operator_export1, libfly), Nothing, (Ptr{Cvoid}, Cint, Float32, Ptr{Float32}, Ptr{Float32}), self.ptr, arg0, newarg1, real, imag)
#    return complex(real[], imag[])
#end

function (self::Field_C)(arg0::Vector{Float64}, arg1=0.0)::ComplexF32
    newarg0 = Float32.(arg0)
    lenarg0 = length(arg0)
    newarg1 = Float32(arg1)
    real = Ref{Float32}()
    imag = Ref{Float32}()
    ccall((:Field_C_operator_export2, libfly), Nothing, (Ptr{Cvoid}, Ptr{Float32}, Cint, Float32, Ptr{Float32}, Ptr{Float32}), self.ptr, newarg0, lenarg0, newarg1, real, imag)
    return complex(real[], imag[])
end

function (self::Field_C)(points::Vector{Vector{Float64}}, w=0.0)::Vector{ComplexF32}
    num_points = length(points)
    if num_points == 0
        return ComplexF32[]
    end

    len = length(points[1])
    points_flat = Float32[]
    for p in points
        append!(points_flat, Float32.(p))
    end

    real_output = zeros(Float32, num_points)
    imag_output = zeros(Float32, num_points)
    ccall((:Field_C_operator_export_list, libfly), Cvoid,
          (Ptr{Cvoid}, Ptr{Float32}, Cint, Cint, Float32, Ptr{Float32}, Ptr{Float32}),
          self.ptr, points_flat, num_points, len, Float32(w), real_output, imag_output)
    return complex.(real_output, imag_output)
end

# COMMENTED OUT - no corresponding export
#function (self::Field_C)(arg0::Int, arg1::Vector{Float64}, arg2=0.0)::ComplexF32
#    newarg1 = Float32.(arg1)
#    lenarg1 = length(arg1)
#    newarg2 = Float32(arg2)
#    real = Ref{Float32}()
#    imag = Ref{Float32}()
#    ccall((:Field_C_operator_export3, libfly), Nothing, (Ptr{Cvoid}, Cint, Ptr{Float32}, Cint, Float32, Ptr{Float32}, Ptr{Float32}), self.ptr, arg0, newarg1, lenarg1, newarg2, real, imag)
#    return complex(real[], imag[])
#end

function destroy!(csf::Field_C)
    ccall((:destroy_Field_C, libfly), Cvoid, (Ptr{Cvoid},), csf.cmf)
end

function Base.finalize(obj::Field_C)
    destroy!(obj)
end

function get_data(obj::Field_C)
    """Get data array reshaped in (w,k) format."""
    ptr = ccall((:Field_C_get_data, libfly), Ptr{Cvoid}, (Ptr{Cvoid},), obj.ptr)
    bd = _basedata_from_ptr(ptr)
    data_array = get_data(bd)

    # Reshape from (nk*nw) to (nw, nk)
    if bd.n_indices == 2
        # Matrix: (dim, dim, nk*nw) -> (dim, dim, nk, nw) -> (dim, dim, nw, nk)
        data_reshaped = reshape(data_array, bd.dim_indices, bd.dim_indices, bd.nk, bd.nw)
        return permutedims(data_reshaped, (1, 2, 4, 3))  # swap k and w axes
    else
        # Scalar: (nk*nw,) -> (nk, nw) -> (nw, nk)
        data_reshaped = reshape(data_array, bd.nk, bd.nw)
        return permutedims(data_reshaped, (2, 1))  # transpose to (nw, nk)
    end
end

mutable struct Field_RM
    ptr::Ptr{Cvoid}
    dimension::Int
    mesh::Vector{Int}
    domain::Matrix{Float32}
    w_points::Vector{Float32}
end

function Field_RM()
    ptr = ccall((:Field_RM_export0, libfly), Ptr{Cvoid}, ())
    obj = Field_RM(ptr, 0, Int[], zeros(Float32, 0, 0), Float32[])
    _load_field_metadata!(obj)
    return obj
end

function Field_RM(filename::String)
    ptr = ccall((:Field_RM_export2, libfly), Ptr{Cvoid}, (Cstring,), filename)
    obj = Field_RM(ptr, 0, Int[], zeros(Float32, 0, 0), Float32[])
    _load_field_metadata!(obj)
    return obj
end

function _load_field_metadata!(obj::Field_RM)
    # Get dimension
    obj.dimension = ccall((:Field_RM_get_dimension, libfly), Cint, (Ptr{Cvoid},), obj.ptr)

    # Get mesh
    mesh_size = ccall((:Field_RM_get_mesh_size, libfly), Cint, (Ptr{Cvoid},), obj.ptr)
    if mesh_size > 0
        obj.mesh = zeros(Int32, mesh_size)
        ccall((:Field_RM_get_mesh, libfly), Cvoid, (Ptr{Cvoid}, Ptr{Cint}), obj.ptr, obj.mesh)
        obj.mesh = Int.(obj.mesh)
    end

    # Get domain
    domain_rows = ccall((:Field_RM_get_domain_rows, libfly), Cint, (Ptr{Cvoid},), obj.ptr)
    domain_cols = ccall((:Field_RM_get_domain_cols, libfly), Cint, (Ptr{Cvoid},), obj.ptr)
    if domain_rows > 0 && domain_cols > 0
        domain_flat = zeros(Float32, domain_rows * domain_cols)
        ccall((:Field_RM_get_domain, libfly), Cvoid, (Ptr{Cvoid}, Ptr{Float32}), obj.ptr, domain_flat)
        obj.domain = reshape(domain_flat, domain_cols, domain_rows)'
    end

    # Get w_points
    w_points_size = ccall((:Field_RM_get_w_points_size, libfly), Cint, (Ptr{Cvoid},), obj.ptr)
    if w_points_size > 0
        obj.w_points = zeros(Float32, w_points_size)
        ccall((:Field_RM_get_w_points, libfly), Cvoid, (Ptr{Cvoid}, Ptr{Float32}), obj.ptr, obj.w_points)
    end
end

function (self::Field_RM)(arg0::Vector{Float64}, arg1=0.0)::Matrix{Float32}
    newarg0 = Float32.(arg0)
    lenarg0 = length(arg0)
    newarg1 = Float32(arg1)

    # Allocate space for matrix results
    max_size = 100
    result = zeros(Float32, max_size * max_size)
    matrix_size = Ref{Cint}(0)

    ccall((:Field_RM_operator_export0, libfly), Cvoid,
          (Ptr{Cvoid}, Ptr{Float32}, Cint, Float32, Ptr{Float32}, Ptr{Cint}),
          self.ptr, newarg0, lenarg0, newarg1, result, matrix_size)

    n = matrix_size[]
    if n == 0
        return Matrix{Float32}(undef, 0, 0)
    end

    # Reshape to matrix
    mat = reshape(result[1:n*n], n, n)

    return mat
end

function (self::Field_RM)(points::Vector{Vector{Float64}}, w=0.0)::Vector{Matrix{Float32}}
    num_points = length(points)
    if num_points == 0
        return Matrix{Float32}[]
    end

    len = length(points[1])
    points_flat = Float32[]
    for p in points
        append!(points_flat, Float32.(p))
    end

    max_size = 100
    output = zeros(Float32, num_points * max_size * max_size)
    matrix_size = Ref{Cint}(0)

    ccall((:Field_RM_operator_export_list, libfly), Cvoid,
          (Ptr{Cvoid}, Ptr{Float32}, Cint, Cint, Float32, Ptr{Float32}, Ptr{Cint}),
          self.ptr, points_flat, num_points, len, Float32(w), output, matrix_size)

    n = matrix_size[]
    if n == 0
        return [Matrix{Float32}(undef, 0, 0) for _ in 1:num_points]
    end

    # Reshape to array of matrices
    results = Vector{Matrix{Float32}}(undef, num_points)
    for i in 1:num_points
        start_idx = (i-1) * n * n + 1
        end_idx = i * n * n
        results[i] = reshape(output[start_idx:end_idx], n, n)
    end

    return results
end

function Base.finalize(obj::Field_RM)
    destroy!(obj)
end

function get_data(obj::Field_RM)
    """Get data array reshaped in (w,k) format."""
    ptr = ccall((:Field_RM_get_data, libfly), Ptr{Cvoid}, (Ptr{Cvoid},), obj.ptr)
    bd = _basedata_from_ptr(ptr)
    data_array = get_data(bd)

    # Reshape from (nk*nw) to (nw, nk)
    if bd.n_indices == 2
        # Matrix: (dim, dim, nk*nw) -> (dim, dim, nk, nw) -> (dim, dim, nw, nk)
        data_reshaped = reshape(data_array, bd.dim_indices, bd.dim_indices, bd.nk, bd.nw)
        return permutedims(data_reshaped, (1, 2, 4, 3))  # swap k and w axes
    else
        # Scalar: (nk*nw,) -> (nk, nw) -> (nw, nk)
        data_reshaped = reshape(data_array, bd.nk, bd.nw)
        return permutedims(data_reshaped, (2, 1))  # transpose to (nw, nk)
    end
end

mutable struct Field_CM
    ptr::Ptr{Cvoid}
    dimension::Int
    mesh::Vector{Int}
    domain::Matrix{Float32}
    w_points::Vector{Float32}
end

function Field_CM()
    ptr = ccall((:Field_CM_export0, libfly), Ptr{Cvoid}, ())
    obj = Field_CM(ptr, 0, Int[], zeros(Float32, 0, 0), Float32[])
    _load_field_metadata!(obj)
    return obj
end

function Field_CM(filename::String)
    ptr = ccall((:Field_CM_export2, libfly), Ptr{Cvoid}, (Cstring,), filename)
    obj = Field_CM(ptr, 0, Int[], zeros(Float32, 0, 0), Float32[])
    _load_field_metadata!(obj)
    return obj
end

function _load_field_metadata!(obj::Field_CM)
    # Get dimension
    obj.dimension = ccall((:Field_CM_get_dimension, libfly), Cint, (Ptr{Cvoid},), obj.ptr)

    # Get mesh
    mesh_size = ccall((:Field_CM_get_mesh_size, libfly), Cint, (Ptr{Cvoid},), obj.ptr)
    if mesh_size > 0
        obj.mesh = zeros(Int32, mesh_size)
        ccall((:Field_CM_get_mesh, libfly), Cvoid, (Ptr{Cvoid}, Ptr{Cint}), obj.ptr, obj.mesh)
        obj.mesh = Int.(obj.mesh)
    end

    # Get domain
    domain_rows = ccall((:Field_CM_get_domain_rows, libfly), Cint, (Ptr{Cvoid},), obj.ptr)
    domain_cols = ccall((:Field_CM_get_domain_cols, libfly), Cint, (Ptr{Cvoid},), obj.ptr)
    if domain_rows > 0 && domain_cols > 0
        domain_flat = zeros(Float32, domain_rows * domain_cols)
        ccall((:Field_CM_get_domain, libfly), Cvoid, (Ptr{Cvoid}, Ptr{Float32}), obj.ptr, domain_flat)
        obj.domain = reshape(domain_flat, domain_cols, domain_rows)'
    end

    # Get w_points
    w_points_size = ccall((:Field_CM_get_w_points_size, libfly), Cint, (Ptr{Cvoid},), obj.ptr)
    if w_points_size > 0
        obj.w_points = zeros(Float32, w_points_size)
        ccall((:Field_CM_get_w_points, libfly), Cvoid, (Ptr{Cvoid}, Ptr{Float32}), obj.ptr, obj.w_points)
    end
end

function (self::Field_CM)(arg0::Vector{Float64}, arg1=0.0)::Matrix{ComplexF32}
    newarg0 = Float32.(arg0)
    lenarg0 = length(arg0)
    newarg1 = Float32(arg1)

    # Allocate space for matrix results
    max_size = 100
    real_result = zeros(Float32, max_size * max_size)
    imag_result = zeros(Float32, max_size * max_size)
    matrix_size = Ref{Cint}(0)

    ccall((:Field_CM_operator_export0, libfly), Cvoid,
          (Ptr{Cvoid}, Ptr{Float32}, Cint, Float32, Ptr{Float32}, Ptr{Float32}, Ptr{Cint}),
          self.ptr, newarg0, lenarg0, newarg1, real_result, imag_result, matrix_size)

    n = matrix_size[]
    if n == 0
        return Matrix{ComplexF32}(undef, 0, 0)
    end

    # Reshape to matrix
    real_mat = reshape(real_result[1:n*n], n, n)
    imag_mat = reshape(imag_result[1:n*n], n, n)

    return complex.(real_mat, imag_mat)
end

function (self::Field_CM)(points::Vector{Vector{Float64}}, w=0.0)::Vector{Matrix{ComplexF32}}
    num_points = length(points)
    if num_points == 0
        return Matrix{ComplexF32}[]
    end

    len = length(points[1])
    points_flat = Float32[]
    for p in points
        append!(points_flat, Float32.(p))
    end

    max_size = 100
    real_output = zeros(Float32, num_points * max_size * max_size)
    imag_output = zeros(Float32, num_points * max_size * max_size)
    matrix_size = Ref{Cint}(0)

    ccall((:Field_CM_operator_export_list, libfly), Cvoid,
          (Ptr{Cvoid}, Ptr{Float32}, Cint, Cint, Float32, Ptr{Float32}, Ptr{Float32}, Ptr{Cint}),
          self.ptr, points_flat, num_points, len, Float32(w), real_output, imag_output, matrix_size)

    n = matrix_size[]
    if n == 0
        return [Matrix{ComplexF32}(undef, 0, 0) for _ in 1:num_points]
    end

    # Reshape to array of matrices
    results = Vector{Matrix{ComplexF32}}(undef, num_points)
    for i in 1:num_points
        start_idx = (i-1) * n * n + 1
        end_idx = i * n * n
        real_mat = reshape(real_output[start_idx:end_idx], n, n)
        imag_mat = reshape(imag_output[start_idx:end_idx], n, n)
        results[i] = complex.(real_mat, imag_mat)
    end

    return results
end

function Base.finalize(obj::Field_CM)
    destroy!(obj)
end

function get_data(obj::Field_CM)
    """Get data array reshaped in (w,k) format."""
    ptr = ccall((:Field_CM_get_data, libfly), Ptr{Cvoid}, (Ptr{Cvoid},), obj.ptr)
    bd = _basedata_from_ptr(ptr)
    data_array = get_data(bd)

    # Reshape from (nk*nw) to (nw, nk)
    if bd.n_indices == 2
        # Matrix: (dim, dim, nk*nw) -> (dim, dim, nk, nw) -> (dim, dim, nw, nk)
        data_reshaped = reshape(data_array, bd.dim_indices, bd.dim_indices, bd.nk, bd.nw)
        return permutedims(data_reshaped, (1, 2, 4, 3))  # swap k and w axes
    else
        # Scalar: (nk*nw,) -> (nk, nw) -> (nw, nk)
        data_reshaped = reshape(data_array, bd.nk, bd.nw)
        return permutedims(data_reshaped, (2, 1))  # transpose to (nw, nk)
    end
end

# Hamiltonian
mutable struct Hamiltonian
    ptr::Ptr{Cvoid}

    function Hamiltonian()
        ptr = ccall((:Hamiltonian_export0, libfly), Ptr{Cvoid}, ())
        if ptr == C_NULL
            error("Failed to initialize Hamiltonian")
        end
        obj = new(ptr)
        finalizer(obj) do x
            ccall((:destroy_Hamiltonian, libfly), Cvoid, (Ptr{Cvoid},), x.ptr)
        end
        return obj
    end
end

function (self::Hamiltonian)(k::Vector{Float64}, w::Float64=0.0)::Matrix{ComplexF32}
    k_arr = Float32.(k)
    k_len = Cint(length(k))
    w_val = Float32(w)

    # Allocate space for matrix results
    max_size = 100
    real_result = zeros(Float32, max_size * max_size)
    imag_result = zeros(Float32, max_size * max_size)
    matrix_size = Ref{Cint}(0)

    ccall((:Hamiltonian_operator_export0, libfly), Cvoid,
          (Ptr{Cvoid}, Ptr{Float32}, Cint, Float32, Ptr{Float32}, Ptr{Float32}, Ptr{Cint}),
          self.ptr, k_arr, k_len, w_val, real_result, imag_result, matrix_size)

    n = matrix_size[]
    if n == 0
        return Matrix{ComplexF32}(undef, 0, 0)
    end

    # Reshape to matrix
    real_mat = reshape(real_result[1:n*n], n, n)
    imag_mat = reshape(imag_result[1:n*n], n, n)

    return complex.(real_mat, imag_mat)
end

# End Functions

function save_field_to_file!(path::String, field)
    ccall((:cmf_save, libfly), Cvoid, (Cstring, Ptr{Cvoid}), path, field.cmf)
end

function save_data!(
    path::String,
    data::AbstractArray,
    mesh::AbstractVector{<:Integer} = Int[],
    domain::AbstractMatrix{<:AbstractFloat} = zeros(Float64, 0, 0),
    w_points::AbstractVector{<:AbstractFloat} = Float64[],
    n_indices::Integer = 1,
    dim_indices::Integer = 1,
)
    is_complex = eltype(data) <: Complex

    # Determine if data is matrix or scalar (ignoring vector for now)
    if dim_indices > 1
        # Matrix data
        # Calculate number of matrices based on data shape
        data_shape = size(data)
        num_matrices = prod(data_shape[1:end-2])  # All dimensions except last two (matrix dimensions)
        mat_dim = dim_indices
        save_data_matrix(path, data, num_matrices, mat_dim, is_complex, mesh, domain, w_points)
    else
        # Scalar data
        save_data_scalar(path, data, is_complex, mesh, domain, w_points)
    end
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

function interleave_complex(A::AbstractArray)
    # Julia uses column-major, C++ uses row-major, so transpose before flattening
    if ndims(A) == 2
        A = permutedims(A, (2, 1))
    elseif ndims(A) == 3
        # For 3D: reverse all dimensions to convert column-major to row-major
        A = permutedims(A, (3, 2, 1))
    end

    out = Vector{Float32}(undef, 2 * length(A))
    @inbounds for i in eachindex(A)
        out[2i - 1] = real(A[i])
        out[2i]     = imag(A[i])
    end
    return out
end

function save_field!(filename::String, values, domain, mesh, w_points = [])
    with_w = length(w_points) > 0 ? true : false
    found_mesh = collect(size(values))
    nbnd = 1
    if with_w && found_mesh[1] != length(w_points) || !with_w && found_mesh[1] != mesh[1]
        nbnd = found_mesh[1]
    end
    domain = reshape(domain, :)
    is_complex = eltype(values) <: Complex
    is_vector = false # not yet enabled vector support
    with_n = nbnd > 1
    if !is_complex
        values_c = reshape(values, :)  # 1D view, memory-compatible
    else
        values_c = interleave_complex(ComplexF32.(values))
    end
    ccall((:field_save_export0, libfly), Cvoid, (Cstring, Ptr{Float32}, Ptr{Cint}, Cint, Cint, Ptr{Float32}, Cint, Bool, Bool, Bool, Bool, Ptr{Float32}), filename, Float32.(domain), Cint.(mesh), length(mesh), nbnd, w_points, length(w_points), is_complex, is_vector, with_w, with_n, values_c)
end

function load_config!(path::String)
    ccall((:load_config_export0, libfly), Cvoid, (Cstring,), path)
end

function flatten_real(data::AbstractArray{<:Real})
    # Julia uses column-major, C++ uses row-major, so transpose before flattening
    if ndims(data) == 2
        return Float32.(vec(permutedims(data, (2, 1))))
    elseif ndims(data) == 3
        # For 3D: reverse all dimensions to convert column-major to row-major
        return Float32.(vec(permutedims(data, (3, 2, 1))))
    else
        return Float32.(vec(data))
    end
end

# Save data functions
function save_data_scalar(filename::String, data::AbstractArray,
                          is_complex::Bool, mesh::Vector{<:Integer}, domain::Matrix{<:Real},
                          w_points::Vector{<:Real}=Float32[])
    # Flatten and interleave data
    if is_complex
        data_interleaved = interleave_complex(data)
    else
        data_interleaved = flatten_real(data)
    end

    total_size = length(vec(data))
    mesh_i32 = Int32.(mesh)
    mesh_size = length(mesh_i32)
    domain_f32 = Float32.(domain)
    domain_rows, domain_cols = size(domain_f32)
    w_points_f32 = Float32.(w_points)
    w_size = length(w_points_f32)

    # Flatten domain
    domain_flat = reshape(domain_f32', :)

    ccall((:save_data_scalar_export0, libfly), Cvoid,
          (Cstring, Ptr{Float32}, Cint, Bool,
           Ptr{Cint}, Cint, Ptr{Float32}, Cint, Cint, Ptr{Float32}, Cint),
          filename, data_interleaved, total_size, is_complex,
          mesh_i32, mesh_size, domain_flat, domain_rows, domain_cols, w_points_f32, w_size)
end

function save_data_vector(filename::String, data::AbstractArray,
                          nk::Integer, vec_len::Integer, is_complex::Bool,
                          mesh::Vector{<:Integer}, domain::Matrix{<:Real},
                          w_points::Vector{<:Real}=Float32[])
    # Flatten and interleave data
    if is_complex
        data_interleaved = interleave_complex(data)
    else
        data_interleaved = flatten_real(data)
    end

    nk_i32 = Int32(nk)
    vec_len_i32 = Int32(vec_len)
    mesh_i32 = Int32.(mesh)
    mesh_size = length(mesh_i32)
    domain_f32 = Float32.(domain)
    domain_rows, domain_cols = size(domain_f32)
    w_points_f32 = Float32.(w_points)
    w_size = length(w_points_f32)

    # Flatten domain
    domain_flat = reshape(domain_f32', :)

    ccall((:save_data_vector_export0, libfly), Cvoid,
          (Cstring, Ptr{Float32}, Cint, Cint, Bool,
           Ptr{Cint}, Cint, Ptr{Float32}, Cint, Cint, Ptr{Float32}, Cint),
          filename, data_interleaved, nk_i32, vec_len_i32, is_complex,
          mesh_i32, mesh_size, domain_flat, domain_rows, domain_cols, w_points_f32, w_size)
end

function save_data_matrix(filename::String, data::AbstractArray,
                          num_matrices::Integer, mat_dim::Integer, is_complex::Bool,
                          mesh::Vector{<:Integer}, domain::Matrix{<:Real},
                          w_points::Vector{<:Real}=Float32[])
    # Flatten and interleave data
    if is_complex
        data_interleaved = interleave_complex(data)
    else
        data_interleaved = flatten_real(data)
    end

    num_matrices_i32 = Int32(num_matrices)
    mat_dim_i32 = Int32(mat_dim)
    mesh_i32 = Int32.(mesh)
    mesh_size = length(mesh_i32)
    domain_f32 = Float32.(domain)
    domain_rows, domain_cols = size(domain_f32)
    w_points_f32 = Float32.(w_points)
    w_size = length(w_points_f32)

    # Flatten domain
    domain_flat = reshape(domain_f32', :)

    ccall((:save_data_matrix_export0, libfly), Cvoid,
          (Cstring, Ptr{Float32}, Cint, Cint, Bool,
           Ptr{Cint}, Cint, Ptr{Float32}, Cint, Cint, Ptr{Float32}, Cint),
          filename, data_interleaved, num_matrices_i32, mat_dim_i32, is_complex,
          mesh_i32, mesh_size, domain_flat, domain_rows, domain_cols, w_points_f32, w_size)
end

# BaseData exports
mutable struct BaseData
    ptr::Ptr{Cvoid}
    is_complex::Bool
    is_vector::Bool
    is_matrix::Bool
    with_k::Bool
    with_w::Bool
    as_mesh::Bool
    n_indices::Int32
    dim_indices::Int32
    dimension::Int32
    nk::Int32
    nw::Int32
    mesh::Vector{Int32}
    domain::Matrix{Float32}
    w_points::Vector{Float32}

    function BaseData(filename::String, ordering::String="k-w")
        # Load from file
        if ordering == "k-w"
            ptr = ccall((:BaseData_load, libfly), Ptr{Cvoid}, (Cstring,), filename)
        else
            ptr = ccall((:BaseData_load_with_ordering, libfly), Ptr{Cvoid},
                       (Cstring, Cstring), filename, ordering)
        end

        if ptr == C_NULL
            error("Failed to load BaseData from file: $filename")
        end

        # Create instance
        obj = new(ptr)

        # Load metadata
        obj.is_complex = Bool(ccall((:BaseData_get_is_complex, libfly), Cint, (Ptr{Cvoid},), ptr))
        obj.is_vector = Bool(ccall((:BaseData_get_is_vector, libfly), Cint, (Ptr{Cvoid},), ptr))
        obj.is_matrix = Bool(ccall((:BaseData_get_is_matrix, libfly), Cint, (Ptr{Cvoid},), ptr))
        obj.with_k = Bool(ccall((:BaseData_get_with_k, libfly), Cint, (Ptr{Cvoid},), ptr))
        obj.with_w = Bool(ccall((:BaseData_get_with_w, libfly), Cint, (Ptr{Cvoid},), ptr))
        obj.as_mesh = Bool(ccall((:BaseData_get_as_mesh, libfly), Cint, (Ptr{Cvoid},), ptr))
        obj.n_indices = ccall((:BaseData_get_n_indices, libfly), Cint, (Ptr{Cvoid},), ptr)
        obj.dim_indices = ccall((:BaseData_get_dim_indices, libfly), Cint, (Ptr{Cvoid},), ptr)
        obj.dimension = ccall((:BaseData_get_dimension, libfly), Cint, (Ptr{Cvoid},), ptr)
        obj.nk = ccall((:BaseData_get_nk, libfly), Cint, (Ptr{Cvoid},), ptr)
        obj.nw = ccall((:BaseData_get_nw, libfly), Cint, (Ptr{Cvoid},), ptr)

        # Load mesh
        mesh_size = ccall((:BaseData_get_mesh_size, libfly), Cint, (Ptr{Cvoid},), ptr)
        if mesh_size > 0
            mesh_buf = Vector{Int32}(undef, mesh_size)
            ccall((:BaseData_get_mesh, libfly), Cvoid, (Ptr{Cvoid}, Ptr{Cint}), ptr, mesh_buf)
            obj.mesh = mesh_buf
        else
            obj.mesh = Int32[]
        end

        # Load domain
        domain_rows = ccall((:BaseData_get_domain_rows, libfly), Cint, (Ptr{Cvoid},), ptr)
        domain_cols = ccall((:BaseData_get_domain_cols, libfly), Cint, (Ptr{Cvoid},), ptr)
        if domain_rows > 0 && domain_cols > 0
            domain_buf = Vector{Float32}(undef, domain_rows * domain_cols)
            ccall((:BaseData_get_domain, libfly), Cvoid, (Ptr{Cvoid}, Ptr{Float32}), ptr, domain_buf)
            obj.domain = reshape(domain_buf, domain_cols, domain_rows)'  # Transpose for Julia column-major
        else
            obj.domain = Matrix{Float32}(undef, 0, 0)
        end

        # Load w_points
        w_size = ccall((:BaseData_get_w_points_size, libfly), Cint, (Ptr{Cvoid},), ptr)
        if w_size > 0
            w_buf = Vector{Float32}(undef, w_size)
            ccall((:BaseData_get_w_points, libfly), Cvoid, (Ptr{Cvoid}, Ptr{Float32}), ptr, w_buf)
            obj.w_points = w_buf
        else
            obj.w_points = Float32[]
        end

        # Register finalizer to cleanup C++ object
        finalizer(obj) do x
            if x.ptr != C_NULL
                ccall((:destroy_BaseData, libfly), Cvoid, (Ptr{Cvoid},), x.ptr)
                x.ptr = C_NULL
            end
        end

        return obj
    end
end

# Helper function to create BaseData from existing pointer (does not manage lifetime)
function _basedata_from_ptr(ptr::Ptr{Cvoid})
    # Directly create object without calling the constructor
    obj = BaseData(ptr, false, false, false, false, false, false, 0, 0, 0, 0, 0,
                   Int32[], Matrix{Float32}(undef, 0, 0), Float32[])

    # Load metadata
    obj.is_complex = Bool(ccall((:BaseData_get_is_complex, libfly), Cint, (Ptr{Cvoid},), ptr))
    obj.is_vector = Bool(ccall((:BaseData_get_is_vector, libfly), Cint, (Ptr{Cvoid},), ptr))
    obj.is_matrix = Bool(ccall((:BaseData_get_is_matrix, libfly), Cint, (Ptr{Cvoid},), ptr))
    obj.with_k = Bool(ccall((:BaseData_get_with_k, libfly), Cint, (Ptr{Cvoid},), ptr))
    obj.with_w = Bool(ccall((:BaseData_get_with_w, libfly), Cint, (Ptr{Cvoid},), ptr))
    obj.as_mesh = Bool(ccall((:BaseData_get_as_mesh, libfly), Cint, (Ptr{Cvoid},), ptr))
    obj.n_indices = ccall((:BaseData_get_n_indices, libfly), Cint, (Ptr{Cvoid},), ptr)
    obj.dim_indices = ccall((:BaseData_get_dim_indices, libfly), Cint, (Ptr{Cvoid},), ptr)
    obj.dimension = ccall((:BaseData_get_dimension, libfly), Cint, (Ptr{Cvoid},), ptr)
    obj.nk = ccall((:BaseData_get_nk, libfly), Cint, (Ptr{Cvoid},), ptr)
    obj.nw = ccall((:BaseData_get_nw, libfly), Cint, (Ptr{Cvoid},), ptr)

    # Load mesh
    mesh_size = ccall((:BaseData_get_mesh_size, libfly), Cint, (Ptr{Cvoid},), ptr)
    if mesh_size > 0
        mesh_buf = Vector{Int32}(undef, mesh_size)
        ccall((:BaseData_get_mesh, libfly), Cvoid, (Ptr{Cvoid}, Ptr{Cint}), ptr, mesh_buf)
        obj.mesh = mesh_buf
    else
        obj.mesh = Int32[]
    end

    # Load domain
    domain_rows = ccall((:BaseData_get_domain_rows, libfly), Cint, (Ptr{Cvoid},), ptr)
    domain_cols = ccall((:BaseData_get_domain_cols, libfly), Cint, (Ptr{Cvoid},), ptr)
    if domain_rows > 0 && domain_cols > 0
        domain_buf = Vector{Float32}(undef, domain_rows * domain_cols)
        ccall((:BaseData_get_domain, libfly), Cvoid, (Ptr{Cvoid}, Ptr{Float32}), ptr, domain_buf)
        obj.domain = reshape(domain_buf, domain_cols, domain_rows)'
    else
        obj.domain = Matrix{Float32}(undef, 0, 0)
    end

    # Load w_points
    w_size = ccall((:BaseData_get_w_points_size, libfly), Cint, (Ptr{Cvoid},), ptr)
    if w_size > 0
        w_buf = Vector{Float32}(undef, w_size)
        ccall((:BaseData_get_w_points, libfly), Cvoid, (Ptr{Cvoid}, Ptr{Float32}), ptr, w_buf)
        obj.w_points = w_buf
    else
        obj.w_points = Float32[]
    end

    # DO NOT register finalizer - pointer is owned by Field object
    return obj
end

function save!(obj::BaseData, filename::String, ordering::String="k-w")
    """Save BaseData to HDF5 file with specified ordering."""
    if ordering == "k-w"
        ccall((:BaseData_save, libfly), Cvoid, (Ptr{Cvoid}, Cstring), obj.ptr, filename)
    else
        ccall((:BaseData_save_with_ordering, libfly), Cvoid,
              (Ptr{Cvoid}, Cstring, Cstring), obj.ptr, filename, ordering)
    end
end

function get_data(obj::BaseData)
    """Extract data as Julia array."""
    if obj.n_indices == 2
        # Matrix data
        total_size = obj.nk * obj.nw * obj.dim_indices * obj.dim_indices
        real_buf = Vector{Float32}(undef, total_size)
        imag_buf = obj.is_complex ? Vector{Float32}(undef, total_size) : real_buf

        ccall((:BaseData_get_data_matrix, libfly), Cvoid,
              (Ptr{Cvoid}, Ptr{Float32}, Ptr{Float32}), obj.ptr, real_buf, imag_buf)

        if obj.is_complex
            data = complex.(real_buf, imag_buf)
        else
            data = real_buf
        end

        # Reshape to (dim_indices, dim_indices, nk*nw) - Julia column-major
        return reshape(data, obj.dim_indices, obj.dim_indices, obj.nk * obj.nw)
    else
        # Scalar data
        total_size = obj.nk * obj.nw
        real_buf = Vector{Float32}(undef, total_size)
        imag_buf = obj.is_complex ? Vector{Float32}(undef, total_size) : real_buf

        ccall((:BaseData_get_data_scalar, libfly), Cvoid,
              (Ptr{Cvoid}, Ptr{Float32}, Ptr{Float32}), obj.ptr, real_buf, imag_buf)

        if obj.is_complex
            return complex.(real_buf, imag_buf)
        else
            return real_buf
        end
    end
end

end # module Imports
