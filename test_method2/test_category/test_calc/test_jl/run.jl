#!/usr/bin/env julia
"""
Simple test script for run_julia_method2.
Reads from stdin and prints confirmation.
"""

println("=== Julia method2 test script ===")
println("Reading from stdin...")

# Read all input
config_data = read(stdin, String)

if !isempty(config_data)
    println("Received $(length(config_data)) bytes of config data")
    println("First 100 chars:")
    println(first(config_data, min(100, length(config_data))))
else
    println("No input received")
end

println("Julia script completed successfully")
exit(0)
