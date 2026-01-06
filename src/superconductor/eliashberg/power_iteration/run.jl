include("eliashberg.jl")
using .Eliashberg

function run()
    max_phi = eliashberg_convsum()
    return max_phi
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end



