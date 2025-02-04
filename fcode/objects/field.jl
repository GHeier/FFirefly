module CondensedMatterField

using Interpolations
using DelimitedFiles
#using ConvexHull
using LinearAlgebra
using Statistics

struct CMF{}
    points::Matrix{Float64}
    data
    interp::Interpolations.GriddedInterpolation
    w_interp::Interpolations.GriddedInterpolation
    inv_domain::Matrix{Float64}
    dimension::Int
    with_w::Bool
    filled::Bool
end

function CMF(points::Matrix{Float64}, data, dimension, with_w::Bool)::CMF
    domain = extract_domain(points, dimension)
    temp_inv = expand_domain(domain, with_w)
    inv_domain = inv(domain)
    temp_inv = expand_domain(inv_domain, with_w)
    transformed_points = Array([temp_inv * points[i,:] for i in 1:size(points,1)])
    interp, w_interp = calculate_interp(transformed_points, data, dimension)
    println("Interpolations calculated")
    println("domain: ", domain)
    println("inv_domain: ", inv_domain)
    return CMF(points, data, interp, w_interp, inv_domain, dimension, with_w, true)
end

function CMF(filename::String)::CMF
    data, header = readdlm(filename, header=true)
    dimension, p1, is_complex, is_vector = get_info(header)
    points = Float64.(data[:,1:dimension + p1])
    data = data[:,dimension + p1 + 1:end]
    data = reconstruct_data(data, is_complex, is_vector)
    return CMF(points, data, dimension, Bool(p1))
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
    if psize > dimension
        coords_list[end] = integer_list
    end
    coords_tuple = tuple(coords_list...)
    interp = interpolate(coords_tuple, data_nd, Gridded(Linear()))

    #Step 4: Create interpolation object for inv_w calculation
    w_tuple = tuple(w_coords)
    data_integer = reshape(integer_list, length(w_coords))
    w_interp = interpolate(w_tuple, data_integer, Gridded(Linear()))

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
    while i < length(points)
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

function (cmf::CMF)(q::Vector{Float64}, w::Float64 = 0.0)
    q = cmf.inv_domain * q
    if cmf.with_w
        q = [q; cmf.w_interp(w)]
    end
    return cmf.interp(q...)
end

function test_2dim_real_scalar_field()
    points = [0.0 0.0; 1.0 0.0; 0.0 1.0; 1.0 1.0]
    data = [1.0, 2.0, 3.0, 4.0]
    cmf = CMF(points, data, 2, false)
    println(cmf([0.5, 0.5], 0.0))
    return true
end

function test_2dim_complex_scalar_field()
    points = [0.0 0.0; 1.0 0.0; 0.0 1.0; 1.0 1.0]
    data = [1.0 + 1.0*im, 2.0 + 2.0*im, 3.0 + 3.0*im, 4.0 + 4.0*im]
    cmf = CMF(points, data, 2, false)
    println(cmf([0.5, 0.5], 0.0))
end

function test_2dim_complex_scalar_field_with_w()
    points = [0.0 0.0 0.0; 1.0 0.0 0.0; 0.0 1.0 0.0; 1.0 1.0 0.0; 0.0 0.0 1.0; 1.0 0.0 1.0; 0.0 1.0 1.0; 1.0 1.0 1.0]
    data = [1.0 + 1.0*im, 2.0 + 2.0*im, 3.0 + 3.0*im, 4.0 + 4.0*im, 5.0 + 5.0*im, 6.0 + 6.0*im, 7.0 + 7.0*im, 8.0 + 8.0*im]
    cmf = CMF(points, data, 2, true)
    println(cmf([0.5, 0.5], 0.0))
end

end # module

