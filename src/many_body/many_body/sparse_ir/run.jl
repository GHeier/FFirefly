using Firefly
cfg = Firefly.config

# Load relevant variables from the configuration
kmesh = cfg.k_mesh



function run():
    # Main function call goes here
    println("Hello, World! This is a Firefly run with k-mesh:", kmesh)

if abspath(PROGRAM_FILE) == @__FILE__ # Runs on file execution
    run()
end



