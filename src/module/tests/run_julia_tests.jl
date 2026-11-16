# Run the Julia binding test suite and exit nonzero if any test fails.
#
# Used by scripts/test-all.sh:
#   julia --project=jlpkg/Firefly src/module/tests/run_julia_tests.jl

include(joinpath(@__DIR__, "test_module_interface.jl"))
using .TestModuleInterface

const M = TestModuleInterface
test_names = sort([n for n in names(M; all=true)
                   if startswith(String(n), "test_") && getfield(M, n) isa Function])

failed = Symbol[]
for n in test_names
    ok = try
        getfield(M, n)()
    catch e
        println("$n threw: ", e)
        false
    end
    ok === true || push!(failed, n)
end

println("\ntest_module_interface.jl: $(length(test_names) - length(failed))/$(length(test_names)) passed")
for n in failed
    println("  FAILED: $n")
end
exit(isempty(failed) ? 0 : 1)
