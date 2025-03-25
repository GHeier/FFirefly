module CMField

using Interpolations, LinearAlgebra, Statistics
using DelimitedFiles, Printf
#using Base: ccall, @ccallable

export CMF, save_CMF, save_arr_as_meshgrid

struct CMF{}
    points::Matrix{Float64}
    data
    interp::Interpolations.GriddedInterpolation
    w_interp::Interpolations.GriddedInterpolation
    inv_domain::Matrix{Float64}
    dimension::Int
    with_w::Bool
    bounds::Vector{Float64}
    filled::Bool
end

function CMF(points::Matrix{Float64}, data, dimension, with_w::Bool)::CMF
    domain = extract_domain(points, dimension)
    temp_inv = expand_domain(domain, with_w)
    inv_domain = inv(domain)
    temp_inv = expand_domain(inv_domain, with_w)
    transformed_points = Array([temp_inv * points[i,:] for i in 1:size(points,1)])
    bounds = get_bounds(transformed_points, dimension, with_w)
    interp, w_interp = calculate_interp(transformed_points, data, dimension)
    return CMF(points, data, interp, w_interp, inv_domain, dimension, with_w, bounds, true)
end

function CMF(filename::String)::CMF
    data, header = readdlm(filename, header=true)
    dimension, p1, is_complex, is_vector = get_info(header)
    points = Float64.(data[:,1:dimension + p1])
    data = data[:,dimension + p1 + 1:end]
    data = reconstruct_data(data, is_complex, is_vector)
    return CMF(points, data, dimension, Bool(p1))
end

function get_bounds(points, dim, with_w)
    points_matrix = reduce(hcat, points)'  # Convert to a 2D matrix
    xmin, xmax, ymin, ymax, zmin, zmax, wmin, wmax = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    xmin = minimum(points_matrix[:,1])
    xmax = maximum(points_matrix[:,1])
    if dim > 1
        ymin = minimum(points_matrix[:,2])
        ymax = maximum(points_matrix[:,2])
    end
    if dim > 2
        zmin = minimum(points_matrix[:,3])
        zmax = maximum(points_matrix[:,3])
    end
    if with_w
        wmin = minimum(points_matrix[:,end])
        wmax = maximum(points_matrix[:,end])
    end
    return [xmin, xmax, ymin, ymax, zmin, zmax, wmin, wmax]
end

function get_info(header)
    dimension = 0
    p1 = 0
    is_complex = false
    is_vector = false
    for i in 1:length(header)
        if header[i] == "x" || header[i] == "y" || header[i] == "z"
            dimension += 1
        elseif header[i] == "w" || header[i] == "Im(w)" 
            p1 = 1
        end
        if (occursin("Im(f)", header[i]) || occursin("Re(f)", header[i]))
            is_complex = true
        end
        if occursin("fx", header[i])
            is_vector = true
        end
    end
    return dimension, p1, is_complex, is_vector
end

function reconstruct_data(data, is_complex, is_vector)
    if is_vector
        if is_complex
            new_data = [data[:,2*i - 1] + data[:,2*i]*im for i in 1:size(data,2) ÷ 2]
        else
            new_data = [data[:,i] for i in 1:size(data,2)]
        end
    else
        if is_complex
            new_data = data[:,1] + data[:,2]*im
        else
            new_data = data[:,1]
        end
    end
    return new_data
end

function calculate_interp(points, data, dimension)
    points_matrix = reduce(hcat, points)'  # Convert to a 2D matrix
    psize = length(points[1])
    coords_list = Vector{Vector{Float64}}(undef, psize)
    for i in 1:psize
        coords_list[i] = unique(points_matrix[:,i])
    end

    # Step 2: Reshape the data into a 3D array
    npts_list = [length(coords) for coords in coords_list]
    npts_tuple = tuple(npts_list...)
    data_nd = reshape(data, npts_tuple)  # Reshape data into a nD array

    # Step 3: Create the interpolation object
    w_coords = coords_list[end]
    integer_list = [1.0*i for i in 1:size(w_coords,1)]
    #if psize > dimension
    #    println("psize > dimension")
    #    coords_list[end] = integer_list
    #end
    coords_tuple = tuple(coords_list...)
    interp = interpolate(coords_tuple, data_nd, Gridded(Linear()))

    #Step 4: Create interpolation object for inv_w calculation
    w_tuple = tuple(w_coords)
    data_integer = reshape(integer_list, length(w_coords))
    w_interp = interpolate(w_tuple, data_integer, Gridded(Linear()))
    println(w_tuple)
    println(data_integer)
    println(w_interp(0.0))

    return interp, w_interp
end

function expand_domain(domain::Matrix{Float64}, with_w::Bool)
    if with_w
        new_domain = Matrix{Float64}(undef, size(domain,1) + 1, size(domain,2) + 1)
        new_domain[1:size(domain,1),1:size(domain,2)] = domain
        new_domain[size(domain,1) + 1,1:size(domain,2)] = zeros(size(domain,2))
        new_domain[1:size(domain,1),size(domain,2) + 1] = zeros(size(domain,1))
        new_domain[size(domain,1) + 1,size(domain,2) + 1] = 1
        return new_domain
    else
        return domain
    end
end

function cross_product(a::Vector{Float64}, b::Vector{Float64})
    if length(a) == 3 && length(b) == 3
        return [a[2]*b[3] - a[3]*b[2], a[3]*b[1] - a[1]*b[3], a[1]*b[2] - a[2]*b[1]]
    elseif length(a) == 2 && length(b) == 2
        return a[1]*b[2] - a[2]*b[1]
    else
        return 0
    end
end

function extract_domain(points::Matrix{Float64}, dimension::Int)
    section = 1
    section_sizes = [1, 1, 1]
    e = dimension
    domain = Matrix{Float64}(undef, e, e)
    
    first = points[1,1:e]
    lattice_vec = round.(points[2,1:e] - first, digits=4)
    jump = 1
    for i in 2:size(points,1)
        lattice_vec = round.(points[i,1:e] - first, digits=4)
        if norm(lattice_vec) > 1e-4
            jump = i - 1
            break
        end
    end
    p = first
    prev_p = lattice_vec

    i = 1
    while i < length(points) - 1
        p = points[i+1,1:e]
        current_vec = round.((p - first) / section_sizes[section], digits=4)
        cross = cross_product(lattice_vec, current_vec)
        if norm(cross) < 1e-4
            section_sizes[section] += 1
        else
            domain[end - section + 1,:] = round.(points[i+1-jump,1:e] - first, digits=4)
            jump = i - 1
            i -= 1
            section += 1
            if section >= dimension
                break
            end
            section_sizes[section] += 1
            lattice_vec = round.(p - first, digits=4)
        end
        prev_p = p
        i += jump
    end

    domain[1,:] = round.(points[end - jump + 1,1:e] - first, digits=4)
    return domain
end

function sanitize_within_bounds(value::Float64, min::Float64, max::Float64, tol::Float64)
    if value <= min && value > min - tol
        return min + tol
    elseif value >= max && value < max + tol
        return max - tol
    else
        return value
    end
end

function (cmf::CMF)(q::Vector{Float64}, w::Float64 = 0.0)
    q = cmf.inv_domain * q
    println(q)
    q = [sanitize_within_bounds(q[i], cmf.bounds[2*i-1], cmf.bounds[2*i], 1e-4) for i in 1:length(q)]
    println(q)
    #if cmf.with_w
    #    q = [q; cmf.w_interp(w)]
    #end
    println(q)
    #q = [q; w]
    println(q)
    return cmf.interp(q...)
end

function get_header(data, dimension, with_w)
    is_complex = !isreal(data[1])
    is_vector = (length(data[1]) > 1 && !is_complex || is_complex && length(data[1]) > 2)
    fheader = "f"
    if is_complex
        fheader = "Re(f) Im(f)"
    end
    if is_vector && !is_complex
        fheader = "fx "
    end
    if is_vector && is_complex
        fheader = "Re(fx) Im(fx) "
    end
    header = "x "
    if dimension > 1
        header *= "y "
        if is_vector && !is_complex
            fheader *= "fy "
        elseif is_vector && is_complex
            fheader *= "Re(fy) Im(fy) "
        end
    end
    if dimension > 2
        header *= "z "
        if is_vector && !is_complex
            fheader *= "fz "
        end
        if is_vector && is_complex
            fheader *= "Re(fz) Im(fz) "
        end
    end
    if with_w
        header *= "w "
    end
    header *= fheader
    return header
end

function save_arr_as_meshgrid(points::Matrix{Float64}, data, filename::String, with_w::Bool = false)
    data_r = reduce(hcat, data)'  # Convert to a 2D matrix
    dimension = size(points,2) - with_w
    is_complex = !isreal(data[1])

    if is_complex
        data_r = hcat(real(data_r), imag(data_r))
    end

    data_matrix = hcat(points, data_r)
    header = get_header(data, dimension, with_w)
    header_as_vector = split(header)

    # Format data as strings with 6 decimal places
    formatted_data = [@sprintf("%.6f", x) for x in data_matrix]  # Flatten to vector
    formatted_data = reshape(formatted_data, size(data_matrix))  # Reshape to match matrix size

    # Compute column widths dynamically
    combined_matrix = vcat(reshape(header_as_vector, 1, :), formatted_data)  # Ensure header is a row
    col_widths = [maximum(length.(combined_matrix[:, i])) for i in 1:size(combined_matrix, 2)]

    open(filename, "w") do io
        # Format and write header
        header_str = join([rpad(header_as_vector[i], col_widths[i]) for i in 1:length(header_as_vector)], "  ")
        println(io, header_str)

        # Write formatted data rows
        for row in eachrow(formatted_data)
            println(io, join([rpad(row[i], col_widths[i]) for i in 1:length(row)], "  "))
        end
    end
end

function save_CMF(cmf::CMF, filename::String)
    save_arr_as_meshgrid(cmf.points, cmf.data, filename, cmf.with_w)
end


# C-callable wrapper to create CMF object
#@ccallable function cmf_create(filename::Cstring)::Ptr{CMF}
#    cmf = CMF(unsafe_string(filename))
#    return pointer_from_objref(cmf)
#end
#
## C-callable wrapper for cmf_eval
#@ccallable function cmf_eval_c(cmf_ptr::Ptr{CMF}, q_ptr::Ptr{Cdouble}, w::Cdouble)::Cdouble
#    cmf = unsafe_pointer_to_objref(cmf_ptr)::CMF
#    q = unsafe_wrap(Array, q_ptr, (size(cmf.inv_domain, 2),))  # Adapt to dimensionality
#    return cmf_eval(cmf, q, w)
#end
#
## C-callable function to free memory (not strictly needed unless managing manually)
#@ccallable function cmf_destroy(cmf_ptr::Ptr{CMF})::Cvoid
#    nothing  # Julia's GC will handle it
#end


end # module

