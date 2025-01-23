cmake_minimum_required(VERSION 3.15)
project(agoge LANGUAGES CXX)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Optionally find FFTW and HDF5 libraries
# (You can remove these if you manually specify them later, or if your code doesn't need them)
find_package(HDF5 COMPONENTS CXX HL REQUIRED)

# Add the 'src' subdirectory which contains our source code
add_subdirectory(src)

# (Optional) If you want tests or other directories, add_subdirectory(tests) here
add_subdirectory(problems)