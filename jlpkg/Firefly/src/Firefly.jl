__precompile__(false)

module Firefly

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

# Read in config values
cfg_file = abspath(@__FILE__)[1:end-28] * "build/bin/input.cfg"
if isfile(cfg_file)
    load_config!(cfg_file)
else
    println("No cfg file found. Run fly.x with input to load in variables")
end

end # module Firefly
