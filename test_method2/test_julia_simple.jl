#!/usr/bin/env julia

# Simple test script to verify Julia can read from stdin
println("=== Julia stdin test ===")
config = read(stdin, String)
println("Received $(length(config)) bytes")
println("First 50 chars: ", first(config, min(50, length(config))))
exit(0)
