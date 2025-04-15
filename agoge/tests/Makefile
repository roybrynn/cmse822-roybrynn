CXX = mpicxx
CXXFLAGS = -std=c++17 -Wall -Wextra -O3 -march=native -DNDEBUG
INCLUDES = -I../include/agoge -I../src

# Detect if we can use OpenMP
ifneq ($(shell $(CXX) -fopenmp -E - < /dev/null > /dev/null 2>&1; echo $$?),0)
  # OpenMP not available
  OMPFLAGS =
else
  # OpenMP available
  OMPFLAGS = -fopenmp
endif

CXXFLAGS += $(OMPFLAGS)

SRC_DIR = ../src
AGOGE_SRCS = $(SRC_DIR)/Field3d.cpp $(SRC_DIR)/ParameterSystem.cpp $(SRC_DIR)/EulerSolver.cpp
OBJS = $(patsubst $(SRC_DIR)/%.cpp,%.o,$(AGOGE_SRCS))

.PHONY: all clean test

all: test_field3d_reorder test_euler_solver

test: all
	./test_field3d_reorder
	./test_euler_solver

%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

test_field3d_reorder: test_field3d_reorder.cpp $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $^ -o $@

test_euler_solver: test_euler_solver.cpp $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $^ -o $@

clean:
	rm -f *.o test_field3d_reorder test_euler_solver
