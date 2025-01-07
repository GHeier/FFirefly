"""
sandbox ❯ g++ -std=c++17 -fPIC -shared -o jmod.so jmodule.cpp \                                 13:22:13 
                  -I/home/g/.julia/artifacts/7a508f56099aa725e5f3cd1623d9a33e3787d163/include \
                  -L/home/g/.julia/artifacts/7a508f56099aa725e5f3cd1623d9a33e3787d163/lib \
                  -I$(julia -e 'print(Sys.BINDIR)')/../include/julia \
                  -L$(julia -e 'print(Sys.BINDIR)')/../lib \
                  -lcxxwrap_julia -ljulia

sandbox ❯ set -gx LD_LIBRARY_PATH $LD_LIBRARY_PATH /home/g/.julia/artifacts/fcdf8ae352855a3922b26b2c54db0af1a787fd1e/lib /home/g/.julia/artifacts/7a508f56099aa725e5f3cd1623d9a33e3787d163/lib


sandbox ❯ julia test.jl

"""
# Load the module and generate the functions
module CppHello
  using CxxWrap
  @wrapmodule(() -> joinpath("/home/g/Research/fcode/fcode","jmod"))

  function __init__()
    @initcxx
  end
end

# Call greet and show the result
@show CppHello.add(1, 2)
