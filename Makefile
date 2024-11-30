# Compiler and flags 
CXXFLAGS = -std=c++20 -O3 -fopenmp -g -fPIC
CFLAGS = -O3 -fopenmp -g
FFLAGS = -O3 -fopenmp -g -L /usr/local/lib -I /usr/local/include -ltetrabz 
LIBS = -llapacke -llapack -lblas -lgfortran -lm -lstdc++

# Source files
C_SRC = $(wildcard src/*.c)
CPP_SRC = $(wildcard src/*.cpp)
F90_SRC = $(wildcard src/*.f90)

SRC_FILES = $(C_SRC) $(CPP_SRC) $(F90_SRC)

C_SUBDIRS = $(shell find src/*/ -name '*.c')
CPP_SUBDIRS = $(shell find src/*/ -name '*.cpp')
F90_SUBDIRS = $(shell find src/*/ -name '*.f90')

SUBDIR_FILES = $(C_SUBDIRS) $(CPP_SUBDIRS) $(F90_SUBDIRS)

# Object files (flattened to build/)
C_OBJS = $(patsubst src/%, build/%, $(notdir $(C_SRC:.c=.o)))
CPP_OBJS = $(patsubst src/%, build/%, $(notdir $(CPP_SRC:.cpp=.o)))
F90_OBJS = $(patsubst src/%, build/%, $(notdir $(F90_SRC:.f90=.o)))

SRC_OBJS = $(C_OBJS) $(CPP_OBJS) $(F90_OBJS)

C_SUBDIR_OBJS = $(patsubst src/%, build/%, $(notdir $(C_SUBDIRS:.c=.o)))
CPP_SUBDIR_OBJS = $(patsubst src/%, build/%, $(notdir $(CPP_SUBDIRS:.cpp=.o)))
F90_SUBDIR_OBJS = $(patsubst src/%, build/%, $(notdir $(F90_SUBDIRS:.f90=.o)))

SUBDIR_OBJS = $(C_SUBDIR_OBJS) $(CPP_SUBDIR_OBJS) $(F90_SUBDIR_OBJS)

# Combine all object files
ALL_OBJS = $(SRC_OBJS) $(SUBDIR_OBJS)

# Default target is the main program
all: clear_log fcode.x

# Clear the compile log
clear_log:
	@> compile.log

# Debug output
debug:
	@echo "SRC_FILES: $(SRC_FILES)"
	@echo "SUBDIR_FILES: $(SUBDIR_FILES)"
	@echo "SUBDIRS: $(SUBDIRS)"
	@echo "SRC_OBJS: $(SRC_OBJS)"
	@echo "SUBDIR_OBJS: $(SUBDIR_OBJS)"
	@echo "ALL_OBJS: $(ALL_OBJS)"

# Specify compilation order for specific Fortran files
fortran_compile_order: 
	@gfortran -c src/config_f.f90 -o build/config_f.o $(FFLAGS) -J build 2>>compile.log || { echo "Compilation failed for gap/config_f.f90"; exit 1; }
	@gfortran -c src/response/response_polarization.f90 -o build/response_polarization.o $(FFLAGS) -I build -J build 2>>compile.log || { echo "Compilation failed for gap/response_polarization.f90"; exit 1; }
	@gfortran -c src/response/response.f90 -o build/response.o $(FFLAGS) -I build -J build 2>>compile.log || { echo "Compilation failed for gap/response.f90"; exit 1; }

# Link all object files into the final executable
fcode.x: $(SRC_OBJS) $(SUBDIR_OBJS)
	@mkdir -p build
	@gcc $(SRC_OBJS) $(SUBDIR_OBJS) -o fcode.x $(FFLAGS) $(LIBS) 2>>compile.log || { echo "Build failed"; exit 1; }

build/%.o: src/%.c
	@mkdir -p build
	@gcc -c $< -o $@ $(CFLAGS) 2>>compile.log || { echo "Compilation failed for $<"; exit 1; }

build/%.o: src/%.cpp
	@mkdir -p build
	@g++ -c $< -o $@ $(CXXFLAGS) 2>>compile.log || { echo "Compilation failed for $<"; exit 1; }

build/%.o: src/%.f90
	@mkdir -p build
	@gfortran -c $< -o $@ $(FFLAGS) -I build -J build 2>>compile.log || { echo "Compilation failed for $<"; exit 1; }

build/%.o: $(C_SRC)
	@mkdir -p build
	@gcc -c $< -o $@ $(CFLAGS) 2>>compile.log || { echo "Compilation failed for $<"; exit 1; }

build/%.o: $(CPP_SRC)
	@mkdir -p build
	@g++ -c $< -o $@ $(CXXFLAGS) 2>>compile.log || { echo "Compilation failed for $<"; exit 1; }

build/%.o: $(F90_SRC)
	@mkdir -p build
	@gfortran -c $< -o $@ $(FFLAGS) -I build -J build 2>>compile.log || { echo "Compilation failed for $<"; exit 1; }


.PHONY: clean
clean:
	rm -rf build/*.o fcode.x build/tests/*.o test.exe build/*.mod
