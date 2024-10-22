# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++20 -O3 -fopenmp -g -llapacke -llapack -lblas

# Directories
SRC_DIR = gap
BUILD_DIR = build
TARGET = main.exe

# Source files and corresponding object files
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRCS))

# The target to build the final executable
$(TARGET): $(OBJS)
	@$(CXX) $(OBJS) -o $(TARGET) $(CXXFLAGS) >/dev/null 2>&1 || { echo "Build failed"; exit 1; }

# Rule to compile .cpp files into .o files in build/
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(BUILD_DIR)
	@$(CXX) -c $< -o $@ $(CXXFLAGS) >/dev/null 2>&1 || { echo "Compilation failed for $<"; exit 1; }

# Clean rule to remove object files and the executable
.PHONY: clean
clean:
	rm -rf $(BUILD_DIR)/*.o $(TARGET)
