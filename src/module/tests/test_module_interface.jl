module TestModuleInterface

using LinearAlgebra

# Add imports path
#project_root = dirname(dirname(dirname(dirname(abspath(@__FILE__)))))
#import_path = joinpath(project_root, "src", "module", "imports")
#push!(LOAD_PATH, import_path)

using Firefly

# Vec class tests
function test_vec_constructor_empty()
    try
        v = Firefly.Vec()
        return true
    catch e
        println("Vec empty constructor error: ", e)
        return false
    end
end

function test_vec_constructor_args()
    try
        v = Firefly.Vec(Float32(1.0), Float32(2.0), Float32(3.0))
        return abs(v.x - 1.0) < 0.001 && abs(v.y - 2.0) < 0.001 && abs(v.z - 3.0) < 0.001
    catch e
        println("Vec args constructor error: ", e)
        return false
    end
end

function test_vec_getters()
    try
        v = Firefly.Vec(Float32(1.5), Float32(2.5), Float32(3.5), Float32(0.0), Float32(0.0), Int32(3), Int32(0))
        return (abs(v.x - 1.5) < 0.001 &&
                abs(v.y - 2.5) < 0.001 &&
                abs(v.z - 3.5) < 0.001 &&
                v.dimension == 3)
    catch e
        println("Vec getters error: ", e)
        return false
    end
end

# Load config test
function test_load_config()
    try
        config_path = "/home/g/Research/FFirefly/build/bin/input.cfg"
        if isfile(config_path)
            Firefly.load_config!(config_path)
            return true
        end
        return false
    catch e
        println("load_config error: ", e)
        return false
    end
end

# Save data tests
function test_save_data_scalar()
    try
        tmpdir = mktempdir()
        filename = joinpath(tmpdir, "test_scalar.h5")
        data = ones(Float32, 10, 10)
        mesh = Int32[10, 10]
        domain = Float32[1.0 0.0; 0.0 1.0]
        Firefly.save_data_scalar(filename, data, false, mesh, domain)
        result = isfile(filename)
        rm(tmpdir, recursive=true)
        return result
    catch e
        println("save_data_scalar error: ", e)
        return false
    end
end

function test_save_data_scalar_complex()
    try
        tmpdir = mktempdir()
        filename = joinpath(tmpdir, "test_scalar_complex.h5")
        data = ones(ComplexF32, 10, 10)
        mesh = Int32[10, 10]
        domain = Float32[1.0 0.0; 0.0 1.0]
        Firefly.save_data_scalar(filename, data, true, mesh, domain)
        result = isfile(filename)
        rm(tmpdir, recursive=true)
        return result
    catch e
        println("save_data_scalar complex error: ", e)
        return false
    end
end

function test_save_data_vector()
    try
        tmpdir = mktempdir()
        filename = joinpath(tmpdir, "test_vector.h5")
        nk = 10
        vec_len = 3
        data = ones(Float32, nk, vec_len)
        mesh = Int32[10]
        domain = Float32[1.0;;]
        Firefly.save_data_vector(filename, data, nk, vec_len, false, mesh, domain)
        result = isfile(filename)
        rm(tmpdir, recursive=true)
        return result
    catch e
        println("save_data_vector error: ", e)
        return false
    end
end

function test_save_data_matrix()
    try
        tmpdir = mktempdir()
        filename = joinpath(tmpdir, "test_matrix.h5")
        mat_dim = 2
        num_matrices = 5
        data = ones(Float32, mat_dim, mat_dim, num_matrices)
        mesh = Int32[5]
        domain = Float32[1.0;;]
        Firefly.save_data_matrix(filename, data, num_matrices, mat_dim, false, mesh, domain)
        result = isfile(filename)
        rm(tmpdir, recursive=true)
        return result
    catch e
        println("save_data_matrix error: ", e)
        return false
    end
end

end  # module
