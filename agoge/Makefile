# Makefile for Agoge
# Portable build system for Mac, Linux, Ubuntu, GitHub Codespaces, etc.
# Designed to be simple and flexible for educational purposes.

# ========================
# 1. Compiler and Flags
# ========================

# Allow overriding the compiler via environment variable, else default
CXX := $(CXX)
ifeq ($(CXX),)
    CXX := g++
endif

# Detect Operating System
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    OS := Mac
    CXX := clang++ # On Mac, prefer clang++
else ifeq ($(UNAME_S),Linux)
    OS := Linux
    # CXX remains as set above
else
    OS := Unknown
    # Keep default compiler
endif

# Prefer g++-10 if available
ifeq ($(shell which g++-10),/usr/bin/g++-10)
    CXX := g++-10
endif

# C++ Standard Flag (default to -std=c++20, allow override)
CXX_STD ?= -std=c++2a

# Compiler Flags
CXXFLAGS := $(CXX_STD) -O3 -Wall -Wextra -MMD -MP -g 

# HDF5 Flags using pkg-config
HDF5_CFLAGS := $(shell pkg-config --cflags hdf5_cpp 2>/dev/null)
HDF5_LIBS := $(shell pkg-config --libs hdf5_cpp 2>/dev/null)

# Fallback if pkg-config fails
ifeq ($(strip $(HDF5_CFLAGS)),)
    HDF5_CFLAGS := -I/usr/local/include
endif
ifeq ($(strip $(HDF5_LIBS)),)
    HDF5_LIBS := -lhdf5_cpp -lhdf5
endif

# Include Paths
INCLUDES := -Iinclude/agoge $(HDF5_CFLAGS)

# Linker Flags
LDFLAGS := $(HDF5_LIBS) 

# ========================
# 2. Directories
# ========================

SRC_DIR := src
PROBLEMS_DIR := problems
BUILD_DIR := build
BUILD_SRC_DIR := $(BUILD_DIR)/src
BUILD_PROB_DIR := $(BUILD_DIR)/problems
INCLUDE_DIR := include/agoge

# ========================
# 3. Sources and Objects
# ========================

# Source Files
SRC_FILES := $(wildcard $(SRC_DIR)/*.cpp)
PROB_FILES := $(wildcard $(PROBLEMS_DIR)/*.cpp)

# Object Files
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_SRC_DIR)/%.o,$(SRC_FILES)) \
             $(patsubst $(PROBLEMS_DIR)/%.cpp,$(BUILD_PROB_DIR)/%.o,$(PROB_FILES))

# Dependency Files
DEP_FILES := $(OBJ_FILES:.o=.d)

# ========================
# 4. Targets
# ========================

# Executable Name
TARGET := agoge_run

# Default Target
all: $(TARGET)

# Link Target
$(TARGET): $(OBJ_FILES)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# Compile Source Files
$(BUILD_SRC_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_SRC_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(BUILD_PROB_DIR)/%.o: $(PROBLEMS_DIR)/%.cpp | $(BUILD_PROB_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Create Build Directories
$(BUILD_SRC_DIR):
	mkdir -p $(BUILD_SRC_DIR)

$(BUILD_PROB_DIR):
	mkdir -p $(BUILD_PROB_DIR)

# Include Dependency Files
-include $(DEP_FILES)

# Clean Target
clean:
	rm -rf $(BUILD_DIR) $(TARGET)

# ========================
# 5. Phony Targets
# ========================

.PHONY: all clean
