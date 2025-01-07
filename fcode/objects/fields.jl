module Fields
using Interpolations
using DelimitedFiles

"""
    create_interpolation(file_path::String; dims::Int=3, field_type::Symbol=:scalar, value_type::DataType=Float64)

Create an interpolation field from a data file.

# Arguments:
- `file_path`: Path to the file containing data.
- `dims`: Number of spatial dimensions (default: 3).
- `field_type`: :scalar or :vector to define the type of field.
- `value_type`: Float64 or ComplexF64 for real/complex data.

# File Format:
- Scalar Real: x, y, z, f
- Scalar Complex: x, y, z, f_re, f_im
- Vector Real: x, y, z, fx, fy, fz
- Vector Complex: x, y, z, fx_re, fx_im, fy_re, fy_im, fz_re, fz_im

# Returns:
- An interpolation object.
"""
function create_interpolation(file_path::String; dims::Int=3, field_type::Symbol=:scalar, value_type::DataType=Float64)
    # Load data from file
    data = readdlm(file_path)
    
    # Extract coordinates and values
    points = eachcol(data[:, 1:dims])  # First `dims` columns are coordinates
    values = data[:, dims+1:end]       # Remaining columns are values
    
    # Ensure enough columns for specified type
    num_value_cols = size(values, 2)
    
    # Handle Scalar Field
    if field_type == :scalar
        if value_type == ComplexF64
            @assert num_value_cols == 2 "Complex scalar fields require exactly 2 columns for real and imaginary parts."
            complex_values = complex.(values[:, 1], values[:, 2])
            values = reshape(complex_values, map(length, [unique(points[d]) for d in 1:dims])...)
        else
            @assert num_value_cols == 1 "Real scalar fields require exactly 1 column for values."
            values = reshape(values[:, 1], map(length, [unique(points[d]) for d in 1:dims])...)
        end
        itp = interpolate(tuple([unique(points[d]) for d in 1:dims]...), values, Gridded(Linear()))
    
    # Handle Vector Field
    elseif field_type == :vector
        num_components = value_type == ComplexF64 ? div(num_value_cols, 2) : num_value_cols
        @assert (value_type == ComplexF64 && num_value_cols % 2 == 0) || (value_type == Float64) "Mismatch in value columns for vector field."
        
        reshaped_values = []
        for i in 1:num_components
            if value_type == ComplexF64
                component_values = complex.(values[:, 2*i-1], values[:, 2*i])
            else
                component_values = values[:, i]
            end
            push!(reshaped_values, reshape(component_values, map(length, [unique(points[d]) for d in 1:dims])...))
        end
        itp = [interpolate(tuple([unique(points[d]) for d in 1:dims]...), v, Gridded(Linear())) for v in reshaped_values]
    
    else
        error("Invalid field_type. Use :scalar or :vector.")
    end
    
    return itp
end
end
