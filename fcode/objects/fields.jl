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
    - Data should be in columns: x, y, z, ... f (for scalar) or x, y, z, ... fx, fy, fz, ...

    # Returns:
    - An interpolation object.
    """
    function create_interpolation(file_path::String; dims::Int=3, field_type::Symbol=:scalar, value_type::DataType=Float64)
        # Load data from file
        data = readdlm(file_path)
        
        # Extract coordinates and values
        points = eachcol(data[:, 1:dims])  # First `dims` columns are coordinates
        values = data[:, dims+1:end]       # Remaining columns are values
        
        # Reshape data into grid
        unique_coords = [unique(points[d]) for d in 1:dims]
        grid = tuple(unique_coords...)
        
        if field_type == :scalar
            values = reshape(values[:, 1], map(length, grid)...)  # Scalar field
            itp = interpolate(grid, values, Gridded(Linear()))
        elseif field_type == :vector
            value_dim = size(values, 2)
            reshaped_values = [reshape(values[:, i], map(length, grid)...) for i in 1:value_dim]
            itp = [interpolate(grid, v, Gridded(Linear())) for v in reshaped_values]
        else
            error("Invalid field_type. Use :scalar or :vector.")
        end
        
        return itp
    end
end
