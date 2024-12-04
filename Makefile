# Compiler and flags 
CXXFLAGS = -std=c++20 -O3 -fopenmp -g -fPIC
CFLAGS = -fopenmp 
FFLAGS = -O3 -fopenmp -L /usr/local/lib -I /usr/local/include -ltetrabz -I src/config/load/
LIBS = -llapacke -llapack -lblas -lgfortran -lm -lstdc++

C_SRCS = $(shell find src/ -name "*.c")
C_TMP = $(notdir $(C_SRCS))
C_OBJS = $(patsubst %.c, build/%.o, $(C_TMP))

CPP_SRCS = $(shell find src/ -name "*.cpp")
CPP_TMP = $(notdir $(CPP_SRCS))
CPP_OBJS = $(patsubst %.cpp, build/%.o, $(CPP_TMP))

F90_SRCS = $(shell find src/ -name "*.f90")
F90_TMP1 = $(notdir $(F90_SRCS))
F90_TMP2 = $(filter-out fortran_config.f90, $(F90_TMP1))
F90_OBJS = $(patsubst %.f90, build/%.o, $(F90_TMP1))

SRC_OBJS = $(C_OBJS) $(CPP_OBJS) $(F90_OBJS)

VPATH = src/config/load $(dir $(C_SRCS)) $(dir $(CPP_SRCS)) $(dir $(F90_SRCS))

# Default target is the main program
all: clear_log fcode.x

# Clear the compile log
clear_log:
	@mkdir -p build
	@> compile.log

# Debug output
debug:
	@echo "VPATH: $(VPATH)"
	@echo "C_SRCS: $(C_SRCS)"
	@echo "C_OBJS: $(C_OBJS)"
	@echo "CPP_SRCS: $(CPP_SRCS)"
	@echo "CPP_OBJS: $(CPP_OBJS)"
	@echo "F90_SRCS: $(F90_SRCS)"
	@echo "F90_OBJS: $(F90_OBJS)"
	@echo "F90_NO_CFG_OBJS: $(F90_NO_CFG_OBJS)"

# Link all object files into the final executable
fcode.x: $(C_OBJS) $(CPP_OBJS) $(F90_OBJS)
	@gcc $(SRC_OBJS) -o fcode.x $(FFLAGS) $(LIBS) 2>>compile.log || { echo "Build failed"; exit 1; }

build/%.o: %.c
	@gcc -c $< -o $@ $(CFLAGS) 2>>compile.log || { echo "Compilation failed for $<"; exit 1; }

build/%.o: %.cpp
	@g++ -c $< -o $@ $(CXXFLAGS) 2>>compile.log || { echo "Compilation failed for $<"; exit 1; }

build/%.o: %.f90
	@gfortran -c $< -o $@ $(FFLAGS) -I build -J build 2>>compile.log || { echo "Compilation failed for $<"; exit 1; }

#build/%.o: $(F90_SRCS)
#	@gfortran -c $< -o $@ $(FFLAGS) -I build -J build 2>>compile.log || { echo "Compilation failed for $<"; exit 1; }

.PHONY: clean
clean:
	@rm -rf build/*.o fcode.x build/tests/*.o test.exe build/*.mod
