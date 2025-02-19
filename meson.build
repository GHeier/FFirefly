project('FrankensteinCode', 'c', 'cpp', 'fortran',
  version: '1.0',
  default_options: ['cpp_std=c++20']
)

cc = meson.get_compiler('c')
cpp = meson.get_compiler('cpp')
fc = meson.get_compiler('fortran')  # Fortran compiler

ccache = find_program('ccache', required: false)
define = ['_GLIBCXX_USE_CXX11_ABI=0']

# Compiler flags
c_args = ['-fopenmp', '-fPIC']
cpp_args = ['-O3', '-fopenmp', '-fPIC']
fortran_args = ['-O3', '-fopenmp', '-L/usr/local/lib', '-I/usr/local/include', '-ltetrabz']

# Dependencies
lapacke = dependency('lapacke', required: false)
lapack = dependency('lapack', required: false)
blas = dependency('blas', required: false)

# Fix for gfortran detection
gfortran = fc.find_library('gfortran', required: false)

mathlib = cc.find_library('m', required: true)
stdcpp = cpp.find_library('stdc++', required: true)
tetrabz = cc.find_library('tetrabz', required: false)
openblas = cc.find_library('openblas', required: false)
julia = cc.find_library('julia', required: false)

# Find Python
python = import('python').find_installation()
python_dep = dependency('python3', required: true)

# Include directories
inc_dirs_obj = include_directories(
  meson.project_source_root(),
  python.get_variable('include_dirs'),
  '/home/g/.julia/juliaup/julia-1.11.3+0.x64.linux.gnu/include/julia'
)

# Link directories
lib_dirs = [
  '/home/g/.julia/juliaup/julia-1.11.3+0.x64.linux.gnu/lib',
  python.get_variable('libdir')
]

# Collect source files dynamically
fs = import('fs')

# Get all source files recursively, excluding "archive/" directories
all_sources = fs.glob(meson.project_source_root(), '**/*.[c|cpp|f90]')
sources = []
foreach file : all_sources
    if 'archive/' not in file
        sources += file
    endif
endforeach

# Executable
executable('fcode.x',
  sources,
  dependencies: [lapacke, lapack, blas, gfortran, mathlib, stdcpp, tetrabz, openblas, julia, python_dep],
  include_directories: inc_dirs_obj,  # Use the correct variable
  link_args: lib_dirs,
  c_args: c_args,
  cpp_args: cpp_args,
  fortran_args: fortran_args,
  install: true
)
