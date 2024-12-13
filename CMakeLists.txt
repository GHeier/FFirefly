# Minimum CMake version required
cmake_minimum_required(VERSION 3.14)

# Project name and version
project(FrankensteinCode VERSION 1.0 LANGUAGES C CXX Fortran)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_definitions(-D_GLIBCXX_USE_CXX11_ABI=0)

# Output directories for binaries and libraries
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Include directories for the project
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/src/config/load)

# Compiler flags
set(CMAKE_CXX_FLAGS "-fopenmp -g -fPIC")
set(CMAKE_C_FLAGS "-fopenmp")
set(CMAKE_Fortran_FLAGS "-fopenmp -L/usr/local/lib -I/usr/local/include -ltetrabz")

# Libraries to link
set(LIBS -llapacke -llapack -lblas -lgfortran -lm -lstdc++ -ltetrabz -lopenblas)

# Find Python3 with required components
find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
include_directories(${Python3_INCLUDE_DIRS})
link_directories(${Python3_LIBRARY_DIRS})

# Compile other shared libraries (if needed)
# This step remains unchanged for other libraries

# Create the Python module "fcode" with specific dependencies
add_library(fcode MODULE
    ${CMAKE_SOURCE_DIR}/src/config/load/pybind_module.cpp
    ${CMAKE_SOURCE_DIR}/src/hamiltonian/potential.cpp
    ${CMAKE_SOURCE_DIR}/src/hamiltonian/band_structure.cpp
    ${CMAKE_SOURCE_DIR}/src/objects/vec.cpp
)

# Link required libraries to fcode
target_link_libraries(fcode PRIVATE ${Python3_LIBRARIES} ${LIBS})

# Set the output name and directory for fcode.so
set_target_properties(fcode PROPERTIES
    PREFIX ""
    SUFFIX ".so"
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/src)  # Place fcode.so in src/

# Install the Python module in the src/ directory
install(TARGETS fcode DESTINATION ${CMAKE_SOURCE_DIR}/src)

# Main executable "fcode.x"
file(GLOB_RECURSE MAIN_SOURCES ${CMAKE_SOURCE_DIR}/src/*.c
                               ${CMAKE_SOURCE_DIR}/src/*.cpp
                               ${CMAKE_SOURCE_DIR}/src/*.f90)
add_executable(fcode.x ${MAIN_SOURCES})

# Set the output directory explicitly for the executable
set_target_properties(fcode.x PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

# Link all shared libraries and Python explicitly to the main executable
target_link_libraries(fcode.x PRIVATE ${Python3_LIBRARIES} ${LIBS})
