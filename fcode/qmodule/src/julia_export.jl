module julia_export

# Automatically find and include all .jl files in the "src" directory
src_dir = joinpath(@__DIR__, "../..")  # Change to your folder name
for file in filter(f -> endswith(f, ".jl"), readdir(src_dir))
    include(joinpath(src_dir, file))
end

end  # module MyProject

using PackageCompiler
create_sysimage([:julia_export], sysimage_path="precompiled_julia.so")
