using Firefly
cfg = Firefly.Config

# Load relevant variables from the configuration
kmesh = cfg.k_mesh

include("linearized_eliashberg.jl")
using .Linearized_Eliashberg


function run()
    # Main function call goes here
    Linearized_Eliashberg.eigenvalue_computation()
end


if abspath(PROGRAM_FILE) == @__FILE__ # Runs on file execution
    run()
end



