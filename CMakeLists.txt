# Minimum CMake version required
cmake_minimum_required(VERSION 3.14)

# Project name and version
project(FrankensteinCode VERSION 1.0 LANGUAGES C CXX Fortran)

# Skip fortran dyndep creation
#set(CMAKE_Fortran_MODULE_DIRECTORY "")

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_definitions(-D_GLIBCXX_USE_CXX11_ABI=0)

set(CMAKE_CXX_COMPILER_LAUNCHER ccache)
set(CMAKE_C_COMPILER_LAUNCHER ccache)

# Output directories for binaries and libraries
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Include directories for the project
include_directories(${CMAKE_SOURCE_DIR})
#include_directories(${CMAKE_SOURCE_DIR}/fcode/config/load)

# Compiler flags
set(CMAKE_CXX_FLAGS "-O3 -fopenmp -fPIC")
set(CMAKE_C_FLAGS "-fopenmp -fPIC")
set(CMAKE_Fortran_FLAGS "-O3 -fopenmp -L/usr/local/lib -I/usr/local/include -ltetrabz")

# Libraries to link
set(LIBS -llapacke -llapack -lblas -lgfortran -lm -lstdc++ -ltetrabz -lopenblas -ljulia)

# Find Python3 to include its headers and libraries
find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
include_directories(${Python3_INCLUDE_DIRS})
link_directories(${Python3_LIBRARY_DIRS})

# Set Julia include and lib paths
set(JULIA_INCLUDE_DIR "/home/g/.julia/juliaup/julia-1.11.3+0.x64.linux.gnu/include/julia")
set(JULIA_LIB_DIR "/home/g/.julia/juliaup/julia-1.11.3+0.x64.linux.gnu/lib")

# Add the Julia include path to the compiler
include_directories(${JULIA_INCLUDE_DIR})

# Link against the Julia library
link_directories(${JULIA_LIB_DIR})

# Main executable "fcode.x"
file(GLOB_RECURSE MAIN_SOURCES ${CMAKE_SOURCE_DIR}/fcode/*.c
                               ${CMAKE_SOURCE_DIR}/fcode/*.cpp
                               ${CMAKE_SOURCE_DIR}/fcode/*.f90)

# Remove files under any 'archive' directory
foreach(file IN LISTS MAIN_SOURCES)
    if(file MATCHES "/archive/" OR file MATCHES "/spline/")
        list(REMOVE_ITEM MAIN_SOURCES ${file})
    endif()
endforeach()

# Remove fmodule.cpp from the source list if it exists
list(REMOVE_ITEM MAIN_SOURCES ${CMAKE_SOURCE_DIR}/fcode/fmodule.cpp)
list(REMOVE_ITEM MAIN_SOURCES ${CMAKE_SOURCE_DIR}/fcode/jmodule.cpp)

add_executable(fcode.x ${MAIN_SOURCES})

# Set the output directory explicitly for the executable
set_target_properties(fcode.x PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

# Link all shared libraries to the main executable
target_link_libraries(fcode.x PRIVATE ${Python3_LIBRARIES} ${LIBS})
