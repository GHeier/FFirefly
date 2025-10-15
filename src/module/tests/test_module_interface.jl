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

# Round-trip tests - save and read back
function test_save_read_scalar_real()
    try
        tmpdir = mktempdir()
        filename = joinpath(tmpdir, "test_roundtrip_scalar_real.h5")

        # Create test data - 10x10 grid with values = x + y
        nk1, nk2 = 10, 10
        data = zeros(Float32, nk1, nk2)
        for i in 1:nk1
            for j in 1:nk2
                data[i, j] = Float32(i + j - 2)  # -2 because Julia is 1-indexed
            end
        end

        mesh = Int32[nk1, nk2]
        domain = Float32[1.0 0.0; 0.0 1.0]

        # Save the data
        Firefly.save_data_scalar(filename, data, false, mesh, domain)

        # Read it back using Field_R
        field = Firefly.Field_R(filename)

        # Verify data at interior points only (avoid boundaries)
        passed = true
        tolerance = 1e-3
        for i in [2, 4, 6, 8]
            for j in [2, 4, 6, 8]
                k_point = Float64[(i-1)/(nk1-1) - 0.5, (j-1)/(nk2-1) - 0.5, 0.0]
                value = field(k_point)
                expected = data[i, j]
                if abs(value - expected) > tolerance
                    println("Mismatch at ($i,$j): got $value, expected $expected")
                    passed = false
                end
            end
        end

        rm(tmpdir, recursive=true)
        return passed
    catch e
        println("save_read_scalar_real error: ", e)
        Base.show_backtrace(stdout, catch_backtrace())
        return false
    end
end

function test_save_read_scalar_complex()
    try
        tmpdir = mktempdir()
        filename = joinpath(tmpdir, "test_roundtrip_scalar_complex.h5")

        # Create test data - 8x8 grid with complex values
        nk1, nk2 = 8, 8
        data = zeros(ComplexF32, nk1, nk2)
        for i in 1:nk1
            for j in 1:nk2
                data[i, j] = ComplexF32(Float32(i-1), Float32(j-1))  # -1 because Julia is 1-indexed
            end
        end

        mesh = Int32[nk1, nk2]
        domain = Float32[1.0 0.0; 0.0 1.0]

        # Save the data
        Firefly.save_data_scalar(filename, data, true, mesh, domain)

        # Read it back using Field_C
        field = Firefly.Field_C(filename)

        # Verify data at interior points only
        passed = true
        tolerance = 1e-3
        for i in [2, 4, 6]
            for j in [2, 4, 6]
                k_point = Float64[(i-1)/(nk1-1) - 0.5, (j-1)/(nk2-1) - 0.5, 0.0]
                value = field(k_point)
                expected = data[i, j]
                if abs(real(value) - real(expected)) > tolerance || abs(imag(value) - imag(expected)) > tolerance
                    println("Mismatch at ($i,$j): got $value, expected $expected")
                    passed = false
                end
            end
        end

        rm(tmpdir, recursive=true)
        return passed
    catch e
        println("save_read_scalar_complex error: ", e)
        Base.show_backtrace(stdout, catch_backtrace())
        return false
    end
end

function test_save_read_with_frequency()
    try
        tmpdir = mktempdir()
        filename = joinpath(tmpdir, "test_roundtrip_freq.h5")

        # Create test data - 5x5 k-grid with 8 frequencies
        nk1, nk2 = 5, 5
        nw = 8
        data = zeros(Float32, nk1, nk2, nw)
        for i in 1:nk1
            for j in 1:nk2
                for w in 1:nw
                    data[i, j, w] = Float32(i + j + w - 3)  # -3 because Julia is 1-indexed
                end
            end
        end

        mesh = Int32[nk1, nk2]  # Only k-space dimensions
        domain = Float32[1.0 0.0; 0.0 1.0]
        w_points = Float32[Float32(w-1) for w in 1:nw]  # -1 to match 0-indexed

        # Save the data
        Firefly.save_data_scalar(filename, data, false, mesh, domain, w_points)

        # Read it back using Field_R
        field = Firefly.Field_R(filename)

        # Verify data at interior points with different frequencies
        passed = true
        tolerance = 1e-3
        for i in [2, 3, 4]
            for j in [2, 3, 4]
                for w in [2, 4, 6]
                    k_point = Float64[(i-1)/(nk1-1) - 0.5, (j-1)/(nk2-1) - 0.5, 0.0]
                    value = field(k_point, Float64(w-1))
                    expected = data[i, j, w]
                    if abs(value - expected) > tolerance
                        println("Mismatch at ($i,$j,w=$w): got $value, expected $expected")
                        passed = false
                    end
                end
            end
        end

        rm(tmpdir, recursive=true)
        return passed
    catch e
        println("save_read_with_frequency error: ", e)
        Base.show_backtrace(stdout, catch_backtrace())
        return false
    end
end

function test_save_data_dispatcher_real()
    """Test save_data!() automatically dispatches for real data"""
    try
        tmpdir = mktempdir()
        filename = joinpath(tmpdir, "test_dispatcher_real.h5")

        # Create simple real data
        nk1, nk2 = 6, 6
        data = zeros(Float32, nk1, nk2)
        for i in 1:nk1
            for j in 1:nk2
                data[i, j] = Float32(i + j - 2)
            end
        end

        mesh = Int32[nk1, nk2]
        domain = Float32[1.0 0.0; 0.0 1.0]

        # Use save_data! (dispatcher) instead of save_data_scalar
        Firefly.save_data!(filename, data, mesh, domain)

        # Read back and verify (interior points only)
        field = Firefly.Field_R(filename)
        passed = true
        tolerance = 1e-3
        for i in [2, 3, 4, 5]
            for j in [2, 3, 4, 5]
                k_point = Float64[(i-1)/(nk1-1) - 0.5, (j-1)/(nk2-1) - 0.5, 0.0]
                value = field(k_point)
                expected = data[i, j]
                if abs(value - expected) > tolerance
                    println("Dispatcher real mismatch at ($i,$j): got $value, expected $expected")
                    passed = false
                end
            end
        end

        rm(tmpdir, recursive=true)
        return passed
    catch e
        println("save_data_dispatcher_real error: ", e)
        Base.show_backtrace(stdout, catch_backtrace())
        return false
    end
end

function test_save_data_dispatcher_complex()
    """Test save_data!() automatically dispatches for complex data"""
    try
        tmpdir = mktempdir()
        filename = joinpath(tmpdir, "test_dispatcher_complex.h5")

        # Create complex data
        nk1, nk2 = 6, 6
        data = zeros(ComplexF32, nk1, nk2)
        for i in 1:nk1
            for j in 1:nk2
                data[i, j] = ComplexF32(Float32(i-1), Float32(j-1))
            end
        end

        mesh = Int32[nk1, nk2]
        domain = Float32[1.0 0.0; 0.0 1.0]

        # Use save_data! (dispatcher) - should auto-detect complex
        Firefly.save_data!(filename, data, mesh, domain)

        # Read back and verify (interior points only)
        field = Firefly.Field_C(filename)
        passed = true
        tolerance = 1e-3
        for i in [2, 3, 4, 5]
            for j in [2, 3, 4, 5]
                k_point = Float64[(i-1)/(nk1-1) - 0.5, (j-1)/(nk2-1) - 0.5, 0.0]
                value = field(k_point)
                expected = data[i, j]
                if abs(real(value) - real(expected)) > tolerance || abs(imag(value) - imag(expected)) > tolerance
                    println("Dispatcher complex mismatch at ($i,$j): got $value, expected $expected")
                    passed = false
                end
            end
        end

        rm(tmpdir, recursive=true)
        return passed
    catch e
        println("save_data_dispatcher_complex error: ", e)
        Base.show_backtrace(stdout, catch_backtrace())
        return false
    end
end

end  # module
