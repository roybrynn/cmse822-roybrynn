# Compiler and flags
CXX := clang++
CXXFLAGS := -std=c++20 -O3 -Iinclude/agoge -I/Library/Developer/CommandLineTools/usr/include -mmacosx-version-min=10.7 -stdlib=libc++ -Wall -Wextra 
LDFLAGS := -L/usr/lib
LDLIBS := -lhdf5_cpp -lhdf5 # Add more libraries if needed

# Directories
SRC_DIR := src
PROBLEMS_DIR := problems
BUILD_DIR := build
TARGET := agoge_run

# Source files
SOURCES := $(wildcard $(SRC_DIR)/*.cpp) $(wildcard $(PROBLEMS_DIR)/*.cpp)

# Object files
OBJECTS := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/src/%.o,$(wildcard $(SRC_DIR)/*.cpp)) \
           $(patsubst $(PROBLEMS_DIR)/%.cpp,$(BUILD_DIR)/problems/%.o,$(wildcard $(PROBLEMS_DIR)/*.cpp))

# Default target
all: $(TARGET)

# Link the target executable
$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

# Compile source files
$(BUILD_DIR)/src/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)/src
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/problems/%.o: $(PROBLEMS_DIR)/%.cpp | $(BUILD_DIR)/problems
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Create build directories
$(BUILD_DIR)/src:
	mkdir -p $(BUILD_DIR)/src

$(BUILD_DIR)/problems:
	mkdir -p $(BUILD_DIR)/problems

# Clean build artifacts
clean:
	rm -rf $(BUILD_DIR) $(TARGET)

.PHONY: all clean