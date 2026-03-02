#!/usr/bin/env julia
"""Test script for centered parameter in save_data functions."""

include("/home/g/Research/FFirefly/src/module/imports/cpp_imports.jl")

function test_centered()
    # Create test data
    mesh = Int32[4, 4]
    domain = Float32[1.0 0.0; 0.0 1.0]
    data = rand(Float32, 4, 4)

    # Test 1: Save with centered=true (default)
    Imports.save_data!("/tmp/test_centered_true_jl.h5", data, mesh, domain; centered=true)
    bd_true = Imports.BaseData("/tmp/test_centered_true_jl.h5")
    println("Test 1 - centered=true: loaded centered = $(bd_true.centered)")
    @assert bd_true.centered == true "Expected centered=true"

    # Test 2: Save with centered=false
    Imports.save_data!("/tmp/test_centered_false_jl.h5", data, mesh, domain; centered=false)
    bd_false = Imports.BaseData("/tmp/test_centered_false_jl.h5")
    println("Test 2 - centered=false: loaded centered = $(bd_false.centered)")
    @assert bd_false.centered == false "Expected centered=false"

    # Test 3: Default value (should be true)
    Imports.save_data!("/tmp/test_centered_default_jl.h5", data, mesh, domain)
    bd_default = Imports.BaseData("/tmp/test_centered_default_jl.h5")
    println("Test 3 - default: loaded centered = $(bd_default.centered)")
    @assert bd_default.centered == true "Expected default centered=true"

    # Cleanup
    for f in ["/tmp/test_centered_true_jl.h5", "/tmp/test_centered_false_jl.h5", "/tmp/test_centered_default_jl.h5"]
        isfile(f) && rm(f)
    end

    println("\nAll Julia tests passed!")
end

test_centered()
