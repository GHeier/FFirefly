# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++20 -O3 -fopenmp -g -llapacke -llapack -lblas

# Directories
SRC_DIR = gap
BUILD_DIR = build
TEST_DIR = test
TEST_BUILD_DIR = build/tests
TARGET = main.exe
TEST_TARGET = test.exe

# Source files and corresponding object files for main program
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRCS))

# Source files and corresponding object files for test program
TEST_SRCS = $(wildcard $(TEST_DIR)/*.cpp)
TEST_OBJS = $(patsubst $(TEST_DIR)/%.cpp,$(TEST_BUILD_DIR)/%.o,$(TEST_SRCS))

# Default target is the main program
all: $(TARGET)

# The target to build the final executable for the main program
$(TARGET): $(OBJS)
	@$(CXX) $(OBJS) -o $(TARGET) $(CXXFLAGS) >/dev/null 2>&1 || { echo "Build failed"; exit 1; }

# The target to build the final executable for the test program
test: $(TEST_TARGET)

$(TEST_TARGET): $(TEST_OBJS)
	@$(CXX) $(TEST_OBJS) -o $(TEST_TARGET) $(CXXFLAGS) >/dev/null 2>&1 || { echo "Test build failed"; exit 1; }

# Rule to compile .cpp files into .o files in build/ for main program
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(BUILD_DIR)
	@$(CXX) -c $< -o $@ $(CXXFLAGS) >/dev/null 2>&1 || { echo "Compilation failed for $<"; exit 1; }

# Rule to compile .cpp files into .o files in build/tests/ for test program
$(TEST_BUILD_DIR)/%.o: $(TEST_DIR)/%.cpp
	@mkdir -p $(TEST_BUILD_DIR)
	@$(CXX) -c $< -o $@ $(CXXFLAGS) >/dev/null 2>&1 || { echo "Test compilation failed for $<"; exit 1; }

# Clean rule to remove object files and executables
.PHONY: clean
clean:
	rm -rf $(BUILD_DIR)/*.o $(TARGET) $(TEST_BUILD_DIR)/*.o $(TEST_TARGET)
