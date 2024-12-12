# Minimum CMake version required
cmake_minimum_required(VERSION 3.14)

# Project name and version
project(FrankensteinCode VERSION 1.0 LANGUAGES C CXX Fortran)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Output directories for binaries and libraries
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Include directories for the project
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/src/config/load)

# Compiler flags
set(CMAKE_CXX_FLAGS "-O3 -fopenmp -g -fPIC")
set(CMAKE_C_FLAGS "-fopenmp")
set(CMAKE_Fortran_FLAGS "-O3 -fopenmp -L/usr/local/lib -I/usr/local/include -ltetrabz")

# Libraries to link
set(LIBS -llapacke -llapack -lblas -lgfortran -lm -lstdc++ -ltetrabz -lopenblas)

# Find Python3 with required components
find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
include_directories(${Python3_INCLUDE_DIRS})
link_directories(${Python3_LIBRARY_DIRS})

# Collect all source files in the project
file(GLOB_RECURSE C_SOURCES ${CMAKE_SOURCE_DIR}/src/*.c)
file(GLOB_RECURSE CXX_SOURCES ${CMAKE_SOURCE_DIR}/src/*.cpp)
file(GLOB_RECURSE FORTRAN_SOURCES ${CMAKE_SOURCE_DIR}/src/*.f90)

# Separate Fortran sources that define modules and those that depend on them
set(FORTRAN_MODULE_SOURCES "")
set(FORTRAN_DEPENDENT_SOURCES "")

foreach(FSOURCE ${FORTRAN_SOURCES})
    file(STRINGS ${FSOURCE} MODCHECK REGEX "module [a-zA-Z0-9_]+")
    if(MODCHECK)
        list(APPEND FORTRAN_MODULE_SOURCES ${FSOURCE})
    else()
        list(APPEND FORTRAN_DEPENDENT_SOURCES ${FSOURCE})
    endif()
endforeach()

# Compile Fortran module sources into a single library
add_library(fortran_modules SHARED ${FORTRAN_MODULE_SOURCES})
set_target_properties(fortran_modules PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(fortran_modules PRIVATE ${LIBS})

# Compile C and C++ sources into individual shared libraries
foreach(CSOURCE ${C_SOURCES})
    get_filename_component(LIB_NAME ${CSOURCE} NAME_WE)
    add_library(${LIB_NAME} SHARED ${CSOURCE})
    target_link_libraries(${LIB_NAME} PRIVATE ${LIBS})
    set_target_properties(${LIB_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)
endforeach()

foreach(CXXSOURCE ${CXX_SOURCES})
    get_filename_component(LIB_NAME ${CXXSOURCE} NAME_WE)
    add_library(${LIB_NAME} SHARED ${CXXSOURCE})
    target_link_libraries(${LIB_NAME} PRIVATE ${LIBS})
    set_target_properties(${LIB_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)
endforeach()

# Compile Fortran dependent sources into individual shared libraries
foreach(FDEP_SOURCE ${FORTRAN_DEPENDENT_SOURCES})
    get_filename_component(LIB_NAME ${FDEP_SOURCE} NAME_WE)
    add_library(${LIB_NAME} SHARED ${FDEP_SOURCE})
    target_link_libraries(${LIB_NAME} PRIVATE fortran_modules ${LIBS})
    set_target_properties(${LIB_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)
endforeach()

# Create the Python module "fcode"
add_library(fcode MODULE
    ${CMAKE_SOURCE_DIR}/src/config/load/pybind_module.cpp
)

# Link all component libraries to the Python module
foreach(CSOURCE ${C_SOURCES})
    get_filename_component(LIB_NAME ${CSOURCE} NAME_WE)
    target_link_libraries(fcode PRIVATE ${LIB_NAME})
endforeach()

foreach(CXXSOURCE ${CXX_SOURCES})
    get_filename_component(LIB_NAME ${CXXSOURCE} NAME_WE)
    target_link_libraries(fcode PRIVATE ${LIB_NAME})
endforeach()

foreach(FDEP_SOURCE ${FORTRAN_DEPENDENT_SOURCES})
    get_filename_component(LIB_NAME ${FDEP_SOURCE} NAME_WE)
    target_link_libraries(fcode PRIVATE ${LIB_NAME})
endforeach()

target_link_libraries(fcode PRIVATE fortran_modules ${Python3_LIBRARIES} ${LIBS})

# Set the output name to ensure it's a Python module
set_target_properties(fcode PROPERTIES
    PREFIX ""
    SUFFIX ".so")

# Install the Python module in the build directory
install(TARGETS fcode DESTINATION ${CMAKE_BINARY_DIR}/python_package)

# Main executable "fcode.x"
file(GLOB_RECURSE MAIN_SOURCES ${CMAKE_SOURCE_DIR}/src/*.c
                               ${CMAKE_SOURCE_DIR}/src/*.cpp
                               ${CMAKE_SOURCE_DIR}/src/*.f90)
add_executable(fcode.x ${MAIN_SOURCES})

# Link all shared libraries and Python explicitly to the main executable
target_link_libraries(fcode.x PRIVATE fortran_modules ${Python3_LIBRARIES} ${LIBS})
