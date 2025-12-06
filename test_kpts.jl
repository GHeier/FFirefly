#!/usr/bin/env julia
# Test k-point ordering between Python and Julia

nk = 5

# Julia approach
kx = range(-0.5, 0.5, length=nk)
ky = range(-0.5, 0.5, length=nk)
kgrid_julia = [[x, y] for x in kx, y in ky]
kgrid_flat_julia = vec(kgrid_julia)

println("Julia k-points (first 10):")
for i in 1:10
    println("  $i: ", kgrid_flat_julia[i])
end
