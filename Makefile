# Compiler and flags
CXXFLAGS = -std=c++20 -O3 -fopenmp -g -fPIC
CFLAGS = -O3 -fopenmp -g
FFLAGS = -O3 -fopenmp -g -L /usr/local/lib -I /usr/local/include -ltetrabz 
LIBS = -llapacke -llapack -lblas -lgfortran -lm -lstdc++

# Source files and corresponding object files for main program
CPP_SRCS = $(wildcard gap/*.cpp)
C_SRCS = $(wildcard gap/*.c)
F90_SRCS = $(wildcard gap/*.f90)
OBJS = $(patsubst gap/%.cpp,build/%.o, $(CPP_SRCS)) \
       $(patsubst gap/%.c,build/%.o, $(filter-out gap/main.c, $(C_SRCS))) \
       $(patsubst gap/%.f90,build/%.o, $(filter-out gap/config_f.f90 gap/response_polarization.f90 gap/response.f90, $(F90_SRCS)))

# All objects including main.o, but only adding main.o once
ALL_OBJS = $(OBJS) build/config_f.o build/response_polarization.o build/response.o build/main.o

# Source files and corresponding object files for test program
TEST_SRCS = $(wildcard tests/*.cpp)
TEST_OBJS = $(patsubst tests/%.cpp,build/tests/%.o,$(TEST_SRCS))

# Default target is the main program
all: clear_log fortran_compile_order fcode.x

# Target to clear compile.log
clear_log:
	@> compile.log

# Specify compilation order for specific Fortran files
fortran_compile_order: 
	@gfortran -c gap/config_f.f90 -o build/config_f.o $(FFLAGS) -J build 2>>compile.log || { echo "Compilation failed for gap/config_f.f90"; exit 1; }
	@gfortran -c gap/response_polarization.f90 -o build/response_polarization.o $(FFLAGS) -I build -J build 2>>compile.log || { echo "Compilation failed for gap/response_polarization.f90"; exit 1; }
	@gfortran -c gap/response.f90 -o build/response.o $(FFLAGS) -I build -J build 2>>compile.log || { echo "Compilation failed for gap/response.f90"; exit 1; }

# Link all object files into the final executable
fcode.x: $(ALL_OBJS)
	@mkdir -p build
	@gcc $(ALL_OBJS) -o fcode.x $(FFLAGS) $(LIBS) 2>>compile.log || { echo "Build failed"; exit 1; }

# The target to build the final executable for the test program
test: test.exe

test.exe: $(TEST_OBJS) $(OBJS) build/config.o
	@g++ $(TEST_OBJS) $(OBJS) build/config.o -o test.exe $(CXXFLAGS) $(LIBS) 2>>compile.log || { echo "Test build failed"; exit 1; }

# Rule to compile main.c separately
build/main.o: gap/main.c 
	@mkdir -p build
	@gcc -c $< -o $@ $(CFLAGS) 2>>compile.log || { echo "Compilation failed for $<"; exit 1; }

# Rule to compile config.c separately (only once)
build/config.o: gap/config.c
	@mkdir -p build
	@gcc -c $< -o $@ $(CFLAGS) 2>>compile.log || { echo "Compilation failed for $<"; exit 1; }

# Rule to compile .cpp files into .o files in build/ for main program
build/%.o: gap/%.cpp
	@mkdir -p build
	@g++ -c $< -o $@ $(CXXFLAGS) 2>>compile.log || { echo "Compilation failed for $<"; exit 1; }

# Rule to compile .c files into .o files in build/
build/%.o: gap/%.c
	@mkdir -p build
	@gcc -c $< -o $@ $(CFLAGS) 2>>compile.log || { echo "Compilation failed for $<"; exit 1; }

# Rule to compile .cpp files into .o files in build/tests/ for test program
build/tests/%.o: tests/%.cpp
	@mkdir -p build/tests
	@g++ -c $< -o $@ $(CXXFLAGS) 2>>compile.log || { echo "Test compilation failed for $<"; exit 1; }

# Clean rule to remove object files, executables, and .mod files
.PHONY: clean
clean:
	rm -rf build/*.o fcode.x build/tests/*.o test.exe build/*.mod
