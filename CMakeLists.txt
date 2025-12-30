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
include_directories(${CMAKE_SOURCE_DIR}/external/h2pack/include)
include_directories(/home/g/Research/H2Lib/Library)

# Include conda environment headers if CONDA_PREFIX is set
if(DEFINED ENV{CONDA_PREFIX})
    include_directories($ENV{CONDA_PREFIX}/include)
    link_directories($ENV{CONDA_PREFIX}/lib)
    message(STATUS "Using conda environment includes: $ENV{CONDA_PREFIX}/include")
endif()

# Also include system headers for packages installed via pacman
include_directories(/usr/include)
link_directories(/usr/lib)

# Compiler flags
set(CMAKE_CXX_FLAGS "-O3 -fopenmp -fPIC -march=native -ffast-math")
set(CMAKE_C_FLAGS "-fopenmp -fPIC")
set(CMAKE_Fortran_FLAGS "-O3 -fopenmp -fPIC")

# H2Pack library paths
link_directories(${CMAKE_SOURCE_DIR}/external/h2pack/lib)
link_directories(/home/g/Research/H2Lib)

# Libraries to link (removed -ljulia, Python linking handled separately below)
set(LIBS -llapacke -llapack -lblas -lgfortran -lm -lstdc++ -lopenblas -lhdf5_cpp -lhdf5 -lfftw3 -lfftw3f -lh2)

# Python and Julia are NO LONGER linked as libraries
# Python/Julia scripts are executed via shell (run_python_method2/run_julia_method2)
# Find Python3 executable for running scripts
if(DEFINED ENV{CONDA_PREFIX})
    set(Python3_ROOT_DIR $ENV{CONDA_PREFIX})
    message(STATUS "Using Python from conda environment: $ENV{CONDA_PREFIX}")
endif()
find_package(Python3 REQUIRED COMPONENTS Interpreter)
message(STATUS "Python3 executable: ${Python3_EXECUTABLE}")
message(STATUS "Python/Julia NOT linked as libraries - using shell execution instead")

# ---------------------------
# BASE LIBRARY SOURCES
# ---------------------------

# BEGIN AUTO-GENERATED BASE LIBRARY SOURCES
# Base library sources
set(BASE_SOURCES
    src/algorithms/scf_interaction.c
    src/algorithms/integration.cpp
    src/algorithms/linear_algebra.cpp
    src/algorithms/fft.cpp
    src/algorithms/electron_number.cpp
    src/algorithms/interpolate.cpp
    src/algorithms/tests/linalg_tests.cpp
    src/algorithms/tests/fft_tests.cpp
    src/algorithms/tests/all.cpp
    src/config/load/c_config.c
    src/config/load/jl_interface_stub.c
    src/config/load/py_interface_stub.c
    src/config/load/cpp_config.cpp
    src/config/load/tests/py_interface_bool_tests.cpp
    src/config/load/tests/py_interface_typed_tests.cpp
    src/config/load/tests/jl_interface_typed_tests.cpp
    src/config/load/tests/all.cpp
    src/hamiltonian/node.cpp
    src/hamiltonian/models/band_structure.cpp
    src/hamiltonian/models/preloaded_hamiltonians.cpp
    src/hamiltonian/models/interaction.cpp
    src/module/exports/temp.cpp
    src/module/exports/cpp_export.cpp
    src/objects/array.cpp
    src/objects/matrix.cpp
    src/objects/hmatrix.cpp
    src/objects/surfaces.cpp
    src/objects/eigenvec.cpp
    src/objects/vec.cpp
    src/objects/CMField/field.cpp
    src/objects/CMField/hamiltonian.cpp
    src/objects/CMField/base_data.cpp
    src/objects/CMField/field_funcs.cpp
    src/objects/CMField/self_energy.cpp
    src/objects/CMField/data_evaluator.cpp
    src/objects/CMField/bands.cpp
    src/objects/CMField/vertex.cpp
    src/objects/CMField/fields.cpp
    src/objects/tests/hmatrix_tests.cpp
    src/objects/tests/array_tests.cpp
    src/objects/tests/surface_tests.cpp
    src/objects/tests/all.cpp
    src/objects/CMData/cmdata.cpp
    src/objects/CMData/tests/all.cpp
    src/objects/CMField/tests/field_wrapper_tests.cpp
    src/objects/CMField/tests/field_tests.cpp
    src/objects/CMField/tests/base_data_tests.cpp
    src/objects/CMField/tests/band_tests.cpp
    src/objects/CMField/tests/hamiltonian_tests.cpp
    src/objects/CMField/tests/matrix_field_tests.cpp
    src/objects/CMField/tests/tensor_field_tests.cpp
    src/objects/CMField/tests/inds_field_tests.cpp
)

# END AUTO-GENERATED BASE LIBRARY SOURCES

# ---------------------------
# FLY.X SOURCES (main + category nodes + non-C++ methods)
# ---------------------------

# Gather all source files for fly.x
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

# Remove base library sources from ALL_SOURCES to avoid duplication
foreach(file IN LISTS BASE_SOURCES)
    list(REMOVE_ITEM ALL_SOURCES ${CMAKE_SOURCE_DIR}/${file})
endforeach()

# Remove C++ method implementation files (they'll be separate executables)
# This removes src/{category}/{calculation}/{method}/*.cpp for C++ methods
foreach(file IN LISTS ALL_SOURCES)
    # Skip files that are in method directories (3+ levels deep in category dirs)
    # This pattern matches: src/{category}/{calculation}/{method}/file.cpp
    # and also: src/{category}/{calculation}/{method}/tests/file.cpp
    if(file MATCHES "src/(hamiltonian|superconductor|many_body)/[^/]+/[^/]+/")
        # Check if it's not a node.cpp file at the category level
        if(NOT file MATCHES "src/(hamiltonian|superconductor|many_body)/node\\.cpp$")
            list(REMOVE_ITEM ALL_SOURCES ${file})
        endif()
    endif()
endforeach()

# Separate Fortran and C/C++ sources for fly.x
set(FORTRAN_SOURCES "")
set(FLY_CPP_SOURCES "")
set(FLY_C_SOURCES "")

foreach(file IN LISTS ALL_SOURCES)
    if(file MATCHES "\\.f90$")
        list(APPEND FORTRAN_SOURCES ${file})
    elseif(file MATCHES "\\.cpp$")
        list(APPEND FLY_CPP_SOURCES ${file})
    elseif(file MATCHES "\\.c$")
        list(APPEND FLY_C_SOURCES ${file})
    endif()
endforeach()

# ---------------------------
# FORTRAN CONFIG LIBRARY (separate to avoid module conflicts)
# ---------------------------

# Create separate Fortran library to avoid "multiple rules generate ffirefly.mod" error
if(FORTRAN_SOURCES)
    add_library(fortran_config SHARED ${FORTRAN_SOURCES})
    set_target_properties(fortran_config PROPERTIES
        OUTPUT_NAME "fortran_config"
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
        Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/fortran_modules"
    )
endif()

# ---------------------------
# BASE SHARED LIBRARY
# ---------------------------

# Create base.so shared library from core components
add_library(base SHARED ${BASE_SOURCES})
target_link_libraries(base PRIVATE ${LIBS})

# Set SONAME and RPATH for runtime linking
set_target_properties(base PROPERTIES
    OUTPUT_NAME "base"
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    SOVERSION 1
    VERSION 1.0
    BUILD_RPATH "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
    INSTALL_RPATH "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
)

# ---------------------------
# FLY.X EXECUTABLE
# ---------------------------

# Create fly.x executable (link Fortran library instead of compiling sources again)
add_executable(fly.x ${FLY_C_SOURCES} ${FLY_CPP_SOURCES})

# Link fly.x against base.so and fortran_config.so
if(FORTRAN_SOURCES)
    target_link_libraries(fly.x PRIVATE base fortran_config ${LIBS})
else()
    target_link_libraries(fly.x PRIVATE base ${LIBS})
endif()

# Set RPATH so fly.x can find base.so at runtime
set_target_properties(fly.x PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
    BUILD_RPATH "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
    INSTALL_RPATH "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
)

# ---------------------------
# LEGACY FLY SHARED LIBRARY (for Python & Julia bindings)
# ---------------------------

# Keep fly shared library for backward compatibility with Python/Julia
add_library(fly SHARED ${FLY_C_SOURCES} ${FLY_CPP_SOURCES})
if(FORTRAN_SOURCES)
    target_link_libraries(fly PRIVATE base fortran_config ${LIBS})
else()
    target_link_libraries(fly PRIVATE base ${LIBS})
endif()

# ---------------------------
# C++ METHOD EXECUTABLES
# ---------------------------

# BEGIN AUTO-GENERATED CPP METHOD EXECUTABLES
# C++ method executables

# hamiltonian/DOS/tetrahedra
add_executable(hamiltonian_DOS_tetrahedra
    src/hamiltonian/DOS/tetrahedra/run.cpp
)
target_link_libraries(hamiltonian_DOS_tetrahedra PRIVATE base)
set_target_properties(hamiltonian_DOS_tetrahedra PROPERTIES
    OUTPUT_NAME "hamiltonian_DOS_tetrahedra.exe"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    BUILD_RPATH "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
    INSTALL_RPATH "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
)

# hamiltonian/FS/tetrahedra
add_executable(hamiltonian_FS_tetrahedra
    src/hamiltonian/FS/tetrahedra/fs.cpp
    src/hamiltonian/FS/tetrahedra/run.cpp
)
target_link_libraries(hamiltonian_FS_tetrahedra PRIVATE base)
set_target_properties(hamiltonian_FS_tetrahedra PROPERTIES
    OUTPUT_NAME "hamiltonian_FS_tetrahedra.exe"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    BUILD_RPATH "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
    INSTALL_RPATH "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
)

# many_body/vertex/from_susceptibility
add_executable(many_body_vertex_from_susceptibility
    src/many_body/vertex/from_susceptibility/run.cpp
)
target_link_libraries(many_body_vertex_from_susceptibility PRIVATE base)
set_target_properties(many_body_vertex_from_susceptibility PROPERTIES
    OUTPUT_NAME "many_body_vertex_from_susceptibility.exe"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    BUILD_RPATH "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
    INSTALL_RPATH "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
)

# superconductor/bcs/lanczos
add_executable(superconductor_bcs_lanczos
    src/superconductor/bcs/lanczos/run.cpp
)
target_link_libraries(superconductor_bcs_lanczos PRIVATE base)
set_target_properties(superconductor_bcs_lanczos PROPERTIES
    OUTPUT_NAME "superconductor_bcs_lanczos.exe"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    BUILD_RPATH "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
    INSTALL_RPATH "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
)

# superconductor/bcs/power_iteration
add_executable(superconductor_bcs_power_iteration
    src/superconductor/bcs/power_iteration/save_data.cpp
    src/superconductor/bcs/power_iteration/frequency_inclusion.cpp
    src/superconductor/bcs/power_iteration/superconductor.cpp
    src/superconductor/bcs/power_iteration/solver.cpp
    src/superconductor/bcs/power_iteration/run.cpp
    src/superconductor/bcs/power_iteration/utilities.cpp
    src/superconductor/bcs/power_iteration/cfg.cpp
    src/superconductor/bcs/power_iteration/matrix_creation.cpp
)
target_link_libraries(superconductor_bcs_power_iteration PRIVATE base)
set_target_properties(superconductor_bcs_power_iteration PROPERTIES
    OUTPUT_NAME "superconductor_bcs_power_iteration.exe"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    BUILD_RPATH "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
    INSTALL_RPATH "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
)

# END AUTO-GENERATED CPP METHOD EXECUTABLES
