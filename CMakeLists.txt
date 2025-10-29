cmake_minimum_required(VERSION 3.14)

project(Ffirefly VERSION 1.0 LANGUAGES C CXX Fortran)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_definitions(-D_GLIBCXX_USE_CXX11_ABI=1)

# Use ccache to speed up compilation
set(CMAKE_CXX_COMPILER_LAUNCHER ccache)
set(CMAKE_C_COMPILER_LAUNCHER ccache)

# Output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Include directories
include_directories(${CMAKE_SOURCE_DIR})

# Compiler flags
set(CMAKE_CXX_FLAGS "-O3 -fopenmp -fPIC")
set(CMAKE_C_FLAGS "-fopenmp -fPIC")
set(CMAKE_Fortran_FLAGS "-O3 -fopenmp -fPIC -L/usr/local/lib -I/usr/local/include -ltetrabz")

# Libraries to link
set(LIBS -llapacke -llapack -lblas -lgfortran -lm -lstdc++ -ltetrabz -lopenblas -ljulia -lhdf5_cpp -lhdf5)

# Find Python3 headers and libraries
# Use conda environment's Python if CONDA_PREFIX is set
if(DEFINED ENV{CONDA_PREFIX})
    set(Python3_ROOT_DIR $ENV{CONDA_PREFIX})
    message(STATUS "Using Python from conda environment: $ENV{CONDA_PREFIX}")
endif()
find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
message(STATUS "Python3 executable: ${Python3_EXECUTABLE}")
message(STATUS "Python3 include dirs: ${Python3_INCLUDE_DIRS}")
message(STATUS "Python3 libraries: ${Python3_LIBRARIES}")
include_directories(${Python3_INCLUDE_DIRS})
link_directories(${Python3_LIBRARY_DIRS})

# Set Julia include and lib paths
# Updated to use Julia 1.12
set(JULIA_INCLUDE_DIR "/usr/include/julia")
set(JULIA_LIB_DIR "/usr/lib")

include_directories(${JULIA_INCLUDE_DIR})
link_directories(${JULIA_LIB_DIR})

# ---------------------------
# COMMON SOURCE FILES
# ---------------------------

# Gather all source files
file(GLOB_RECURSE ALL_SOURCES 
    ${CMAKE_SOURCE_DIR}/src/*.c
    ${CMAKE_SOURCE_DIR}/src/*.cpp
    ${CMAKE_SOURCE_DIR}/src/*.f90
)

# Remove files under any 'archive' directory
foreach(file IN LISTS ALL_SOURCES)
    if(file MATCHES "/archive/")
        list(REMOVE_ITEM ALL_SOURCES ${file})
    endif()
endforeach()

# Separate Fortran and C/C++ sources
set(FORTRAN_SOURCES "")
set(CPP_SOURCES "")
set(C_SOURCES "")

foreach(file IN LISTS ALL_SOURCES)
    if(file MATCHES "\\.f90$")
        list(APPEND FORTRAN_SOURCES ${file})
      elseif(file MATCHES "\\.cpp$")
        list(APPEND CPP_SOURCES ${file})
      elseif(file MATCHES "\\.c$")
        list(APPEND C_SOURCES ${file})
    endif()
endforeach()

# ---------------------------
# SHARED OBJECT FILES (.o)
# ---------------------------

# Generate object files separately for C/C++ and Fortran
add_library(objlib_c OBJECT ${C_SOURCES})
add_library(objlib_cpp OBJECT ${CPP_SOURCES})
add_library(objlib_fortran OBJECT ${FORTRAN_SOURCES})

# ---------------------------
# SHARED LIBRARY (for Python & Julia) [NO FORTRAN]
# ---------------------------
add_library(fly SHARED $<TARGET_OBJECTS:objlib_c> $<TARGET_OBJECTS:objlib_cpp> $<TARGET_OBJECTS:objlib_fortran>)

# ---------------------------
# EXECUTABLE (Standalone) [WITH FORTRAN]
# ---------------------------
add_executable(fly.x $<TARGET_OBJECTS:objlib_c> $<TARGET_OBJECTS:objlib_cpp> $<TARGET_OBJECTS:objlib_fortran>)

# ---------------------------
# LINKING
# ---------------------------

# Set output directory explicitly
set_target_properties(fly.x PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

# Link necessary libraries (but NOT the shared library!)
target_link_libraries(fly.x PRIVATE ${LIBS} ${Python3_LIBRARIES})
target_link_libraries(fly PRIVATE ${Python3_LIBRARIES} ${LIBS})
