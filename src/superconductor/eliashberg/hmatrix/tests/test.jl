include("../run.jl")
# Set configuration variables for test run
debug = true
mu_from_n = false



function test()
    # Main function call goes here
    result, expected = run()
    println("Expected: ", expected)
    println("Result: ", result)
    return abs(result - expected) < 1e-2  # Example test condition
end


if abspath(PROGRAM_FILE) == @__FILE__ # Runs on file execution
    pass = test()
    if pass
        println("Test passed!")
    else
        println("Test failed.")
    end
end



