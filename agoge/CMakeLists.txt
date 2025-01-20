cmake_minimum_required(VERSION 3.15)
project(agoge LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# If you want to build tests by default
option(AGOGE_BUILD_TESTS "Build agoge unit tests" ON)

# If you use FFTW or HDF5, you can find them here (and add more flags)
find_package(FFTW3 REQUIRED)
find_package(HDF5 COMPONENTS CXX HL REQUIRED)

# Add include directory to be visible project-wide
include_directories(${PROJECT_SOURCE_DIR}/include)

# If you keep your libraries in subdirectories, do:
add_subdirectory(src)

# Optionally, if you keep tests in a separate directory:
if(AGOGE_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
