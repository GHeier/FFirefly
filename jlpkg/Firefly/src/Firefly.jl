__precompile__(false)

module Firefly

export printv

include("./src/config/load/julia_config.jl")

# Handle importing CPP objects and functions
include("./src/module/imports/cpp_imports.jl")
using .Imports
for name in names(Imports, all=false)
    @eval const $(name) = Imports.$(name)
end
for name in names(Imports, all=false)
    @eval export $(name)
end

# Read in config values: $FIREFLY_CONFIG, else build/bin/input.cfg (written by fly.x)
cfg_file = get(ENV, "FIREFLY_CONFIG", "")
if isempty(cfg_file)
    cfg_file = joinpath(dirname(dirname(dirname(@__DIR__))), "build", "bin", "input.cfg")
end
if isfile(cfg_file)
    load_config!(cfg_file)
else
    println("No cfg file found. Run fly.x or set FIREFLY_CONFIG to load in variables")
end

function printv(s)
    if Config.verbosity == "high"
        println(s)
    end
end

end # module Firefly
