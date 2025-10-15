module TestTypedInterface

# Integer tests
function return_int_42()
    return Int32(42)
end

function return_int_negative()
    return Int32(-100)
end

function return_int_zero()
    return Int32(0)
end

# Float tests
function return_float_pi()
    return Float32(3.14159)
end

function return_float_negative()
    return Float32(-2.71828)
end

function return_float_zero()
    return Float32(0.0)
end

# Double tests
function return_double_large()
    return Float64(1.23456789012345)
end

function return_double_small()
    return Float64(0.00000123456)
end

# String tests
function return_string_hello()
    return "Hello, World!"
end

function return_string_empty()
    return ""
end

function return_string_special()
    return "Test!@#\$%^&*()"
end

# Computed values
function compute_sum()
    return Int32(10 + 20 + 30)
end

function compute_product()
    return Float32(2.5 * 4.0)
end

function compute_message()
    name = "FFirefly"
    version = "1.0"
    return "$name v$version"
end

# Boolean test
function return_bool_true()
    return true
end

function return_bool_false()
    return false
end

end  # module

# Test the module
if abspath(PROGRAM_FILE) == @__FILE__
    using .TestTypedInterface

    println("return_int_42(): ", TestTypedInterface.return_int_42())
    println("return_int_negative(): ", TestTypedInterface.return_int_negative())
    println("return_int_zero(): ", TestTypedInterface.return_int_zero())
    println("return_float_pi(): ", TestTypedInterface.return_float_pi())
    println("return_float_negative(): ", TestTypedInterface.return_float_negative())
    println("return_float_zero(): ", TestTypedInterface.return_float_zero())
    println("return_double_large(): ", TestTypedInterface.return_double_large())
    println("return_double_small(): ", TestTypedInterface.return_double_small())
    println("return_string_hello(): ", TestTypedInterface.return_string_hello())
    println("return_string_empty(): '", TestTypedInterface.return_string_empty(), "'")
    println("return_string_special(): ", TestTypedInterface.return_string_special())
    println("compute_sum(): ", TestTypedInterface.compute_sum())
    println("compute_product(): ", TestTypedInterface.compute_product())
    println("compute_message(): ", TestTypedInterface.compute_message())
    println("return_bool_true(): ", TestTypedInterface.return_bool_true())
    println("return_bool_false(): ", TestTypedInterface.return_bool_false())
end
