# Compiler and flags
CXXFLAGS = -std=c++20 -O3 -fopenmp -g
CFLAGS = -O3 -fopenmp -g
FFLAGS = -O3 -fopenmp -g -L /usr/local/lib -I /usr/local/include -ltetrabz
LIBS = -llapacke -llapack -lblas -lgfortran -lm -lstdc++

# Source files and corresponding object files for main program
CPP_SRCS = $(wildcard gap/*.cpp)
C_SRCS = $(wildcard gap/*.c)
F90_SRCS = $(wildcard gap/*.f90)
OBJS = $(patsubst gap/%.cpp,build/%.o, $(CPP_SRCS)) \
       $(patsubst gap/%.c,build/%.o,$(filter-out gap/main.c, $(C_SRCS))) \
       $(patsubst gap/%.f90,build/%.o, $(F90_SRCS))
ALL_OBJS = $(wildcard build/*.o)

# Source files and corresponding object files for test program
TEST_SRCS = $(wildcard tests/*.cpp)
TEST_OBJS = $(patsubst tests/%.cpp,build/tests/%.o,$(TEST_SRCS))

# Default target is the main program
all: main.exe

# Link all object files into the final executable
main.exe: $(ALL_OBJS)
	@gcc $(ALL_OBJS) -o main.exe $(FFLAGS) $(LIBS) 2>build.log || { echo "Build failed"; exit 1; }

# The target to build the final executable for the test program
test: test.exe

test.exe: $(TEST_OBJS) $(OBJS) build/config.o
	@g++ $(TEST_OBJS) $(OBJS) build/config.o -o test.exe $(CXXFLAGS) $(LIBS) 2>build.log || { echo "Test build failed"; exit 1; }

# Rule to compile main.c separately
build/main.o: gap/main.c 
	@mkdir -p build
	@gcc -c $< -o $@ $(CFLAGS) 2>>build.log || { echo "Compilation failed for $<"; exit 1; }

# Rule to compile config.c separately (only once)
build/config.o: gap/config.c
	@mkdir -p build
	@gcc -c $< -o $@ $(CFLAGS) 2>>build.log || { echo "Compilation failed for $<"; exit 1; }

# Rule to compile .cpp files into .o files in build/ for main program
build/%.o: gap/%.cpp
	@mkdir -p build
	@g++ -c $< -o $@ $(CXXFLAGS) 2>>build.log || { echo "Compilation failed for $<"; exit 1; }

# Rule to compile .c files into .o files in build/
build/%.o: gap/%.c
	@mkdir -p build
	@gcc -c $< -o $@ $(CFLAGS) 2>>build.log || { echo "Compilation failed for $<"; exit 1; }

# Rule to compile .f90 files into .o files in build/
build/%.o: gap/%.f90
	@mkdir -p build
	@gfortran -c $< -o $@ $(FFLAGS) 2>>build.log || { echo "Compilation failed for $<"; exit 1; }

# Rule to compile .cpp files into .o files in build/tests/ for test program
build/tests/%.o: tests/%.cpp
	@mkdir -p build/tests
	@g++ -c $< -o $@ $(CXXFLAGS) 2>>build.log || { echo "Test compilation failed for $<"; exit 1; }

# Clean rule to remove object files and executables
.PHONY: clean
clean:
	rm -rf build/*.o main.exe build/tests/*.o test.exe
