module CondensedMatterField

using Interpolations
using DelimitedFiles
#using ConvexHull
using LinearAlgebra
using Statistics

struct CMF{}
    points::Matrix{Float64}
    data::Vector{Float64}
    interp::Interpolations.GriddedInterpolation
    inv_domain::Matrix{Float64}
    dimension::Int
    filled::Bool
end

function CMF(points::Matrix{Float64}, data::Vector{Float64})::CMF
    domain = calculate_domain(points)
    inv_domain = inv(domain)
    transformed_points = Array([inv_domain * points[i,:] for i in 1:size(points,1)])
    interp = calculate_interp(transformed_points, data)
    empty, dimension = size(points)
    return CMF(points, data, interp, inv_domain, dimension, true)
end

function CMF(filename::String)::CMF
    data = readdlm(filename)
    points = data[:,1:end-1]
    data = data[:,end]
    return CMF(points, data)
end

function calculate_interp(points, data)
    points_matrix = reduce(hcat, points)'  # Convert to a 2D matrix
    # Step 1: Extract unique coordinates
    x_coords = unique(points_matrix[:, 1])  # Unique x-coordinates
    y_coords = unique(points_matrix[:, 2])  # Unique y-coordinates
    z_coords = unique(points_matrix[:, 3])  # Unique z-coordinates

    # Step 2: Reshape the data into a 3D array
    nx = length(x_coords)
    ny = length(y_coords)
    nz = length(z_coords)
    data_3d = reshape(data, (nx, ny, nz))  # Reshape data into a 3D array

    # Step 3: Create the interpolation object
    interp = interpolate((x_coords, y_coords, z_coords), data_3d, Gridded(Linear()))

    return interp
end

function calculate_domain(points::Matrix{Float64})
    section = 1
    section_sizes = [1, 1, 1]
    empty, dimension = size(points)
    domain = Matrix{Float64}(undef, dimension, dimension)
    
    first = points[1,:]
    lattice_vec = round.(points[2,:] - first, digits=4)
    jump = 1
    p = first
    prev_p = lattice_vec

    i = 2
    while i <= length(points)
        ind = i
        #if section > 1
        #    ind += 1
        #end
        p = points[ind,:]
        current_vec = round.((p - first) / section_sizes[section], digits=4)
        if norm(current_vec - lattice_vec) < 1e-4
            section_sizes[section] += 1
        else
            domain[section,:] = round.(prev_p - first, digits=4)
            jump *= section_sizes[section]
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

    domain[dimension,:] = round.(points[end - jump + 1,:] - first, digits=4)
    return domain
end

function (cmf::CMF)(q::Vector{Float64})
    q = cmf.inv_domain * q
    return cmf.interp(q...)
end

end # module

using .CondensedMatterField
chi = CondensedMatterField.CMF("/home/g/Research/fcode/chi_mesh_static.dat")
println(chi([0.1, 0.1, 0.1]))
