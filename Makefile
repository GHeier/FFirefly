# Compiler and flags
CXX = g++
CC = gcc
FC = gfortran
CXXFLAGS = -std=c++20 -O3 -fopenmp -g
CFLAGS = -O3 -fopenmp -g
FFLAGS = -O3 -fopenmp -g -L /usr/local/lib -I /usr/local/include -ltetrabz
LIBS = -llapacke -llapack -lblas -lgfortran -lm

# Directories
SRC_DIR = gap
BUILD_DIR = build
TEST_DIR = tests
TEST_BUILD_DIR = build/tests
TARGET = main.exe
TEST_TARGET = test.exe

# Source files and corresponding object files for main program
CPP_SRCS = $(wildcard $(SRC_DIR)/*.cpp)
C_SRCS = $(wildcard $(SRC_DIR)/*.c)
F90_SRCS = $(wildcard $(SRC_DIR)/*.f90)
OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(filter-out $(SRC_DIR)/main.cpp, $(CPP_SRCS))) \
       $(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(filter-out $(SRC_DIR)/main.c $(SRC_DIR)/config.c, $(C_SRCS))) \
       $(patsubst $(SRC_DIR)/%.f90,$(BUILD_DIR)/%.o,$(F90_SRCS))

# Object file for main program
MAIN_OBJ = $(BUILD_DIR)/main.o

# Object file for config (compiled only once)
CONFIG_OBJ = $(BUILD_DIR)/config.o

# Source files and corresponding object files for test program
TEST_SRCS = $(wildcard $(TEST_DIR)/*.cpp)
TEST_OBJS = $(patsubst $(TEST_DIR)/%.cpp,$(TEST_BUILD_DIR)/%.o,$(TEST_SRCS))

# Default target is the main program
all: $(TARGET)

# The target to build the final executable for the main program
$(TARGET): $(OBJS) $(MAIN_OBJ) $(CONFIG_OBJ)
	@$(CC) $(OBJS) $(MAIN_OBJ) $(CONFIG_OBJ) -o $(TARGET) $(CFLAGS) $(LIBS) 2>build.log || { echo "Build failed"; exit 1; }

# The target to build the final executable for the test program
test: $(TEST_TARGET)

$(TEST_TARGET): $(TEST_OBJS) $(OBJS) $(CONFIG_OBJ)
	@$(CXX) $(TEST_OBJS) $(OBJS) $(CONFIG_OBJ) -o $(TEST_TARGET) $(CXXFLAGS) $(LIBS) 2>build.log || { echo "Test build failed"; exit 1; }

# Rule to compile main.c separately
$(BUILD_DIR)/main.o: $(SRC_DIR)/main.c config.h
	@mkdir -p $(BUILD_DIR)
	@$(CC) -c $< -o $@ $(CFLAGS) 2>>build.log || { echo "Compilation failed for $<"; exit 1; }

# Rule to compile config.c separately (only once)
$(BUILD_DIR)/config.o: $(SRC_DIR)/config.c config.h
	@mkdir -p $(BUILD_DIR)
	@$(CC) -c $< -o $@ $(CFLAGS) 2>>build.log || { echo "Compilation failed for $<"; exit 1; }

# Rule to compile .cpp files into .o files in build/ for main program
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(BUILD_DIR)
	@$(CXX) -c $< -o $@ $(CXXFLAGS) 2>>build.log || { echo "Compilation failed for $<"; exit 1; }

# Rule to compile .c files into .o files in build/
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(BUILD_DIR)
	@$(CC) -c $< -o $@ $(CFLAGS) 2>>build.log || { echo "Compilation failed for $<"; exit 1; }

# Rule to compile .f90 files into .o files in build/
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90
	@mkdir -p $(BUILD_DIR)
	@$(FC) -c $< -o $@ $(FFLAGS) 2>>build.log || { echo "Compilation failed for $<"; exit 1; }

# Rule to compile .cpp files into .o files in build/tests/ for test program
$(TEST_BUILD_DIR)/%.o: $(TEST_DIR)/%.cpp
	@mkdir -p $(TEST_BUILD_DIR)
	@$(CXX) -c $< -o $@ $(CXXFLAGS) 2>>build.log || { echo "Test compilation failed for $<"; exit 1; }

# Clean rule to remove object files and executables
.PHONY: clean
clean:
	rm -rf $(BUILD_DIR)/*.o $(TARGET) $(TEST_BUILD_DIR)/*.o $(TEST_TARGET)
