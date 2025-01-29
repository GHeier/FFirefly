#"""
#sandbox ❯ g++ -std=c++17 -fPIC -shared -o jmod.so jmodule.cpp objects/vec.cpp hamiltonian/band_structure.cpp config/load/cpp_config.cpp config/load/c_config.c -I/home/g/.julia/artifacts/7a508f56099aa725e5f3cd1623d9a33e3787d163/include -L/home/g/.julia/artifacts/7a508f56099aa725e5f3cd1623d9a33e3787d163/lib -I$(julia -e 'print(Sys.BINDIR)')/../include/julia -L$(julia -e 'print(Sys.BINDIR)')/../lib -lcxxwrap_julia -ljulia
#
#sandbox ❯ set -gx LD_LIBRARY_PATH $LD_LIBRARY_PATH /home/g/.julia/artifacts/fcdf8ae352855a3922b26b2c54db0af1a787fd1e/lib /home/g/.julia/artifacts/7a508f56099aa725e5f3cd1623d9a33e3787d163/lib
#
#
#sandbox ❯ julia test.jl
#
#"""
module fcode

using CxxWrap

# Wrap C++ shared object
@wrapmodule(() -> joinpath("/home/g/Research/fcode/fcode", "jmod"))
# Initialize C++ environment
function __init__()
    @initcxx
    load_c_config("/home/g/Research/fcode/build/bin/input.cfg")  
end

end # module

