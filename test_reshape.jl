#!/usr/bin/env julia

nk = 3

# Create flat array
vals = collect(0:8)
println("Flat values: ", vals)

# Julia reshape (column-major)
reshaped = reshape(vals, nk, nk)
println("\nJulia reshape (column-major):")
println(reshaped)
println("Element [2,2]: ", reshaped[2,2], " (should be 4)")
