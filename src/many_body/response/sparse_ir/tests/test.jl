include("../run.jl")
# Set configuration variables for test run
kmesh = [4, 4, 4]  # Example k-mesh values



function test()
    # Main function call goes here
    println("Welcome to Testing! This is a Firefly run with k-mesh:", kmesh)
    result = run()
    return abs(result - 3.14) < 1e-6  # Example test condition
end


if abspath(PROGRAM_FILE) == @__FILE__ # Runs on file execution
    run()
end





