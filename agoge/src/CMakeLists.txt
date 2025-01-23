# Go into the 'src' directory for building the core library and the executable
file(GLOB AGOGE_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
    # Alternatively list them explicitly:
    # Field3D.cpp
    # EulerSolver.cpp
    # GravitySolver.cpp
    # HDF5_IO.cpp
    # (main.cpp is built separately for the executable below)
)

# Remove 'main.cpp' from the library if it is included in the wildcard
list(REMOVE_ITEM AGOGE_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp")

# Create a library from the solver source files
add_library(agoge_lib ${AGOGE_SOURCES})

# Add link dependencies
# We'll assume find_package(HDF5 ...) created HDF5::HDF5, HDF5::HDF5_CPP, etc.
target_link_libraries(agoge_lib
    PUBLIC
        HDF5::HDF5
        HDF5::HDF5_CPP
        HDF5::HDF5_HL
)

# Make sure the library can see the headers in this directory
target_include_directories(agoge_lib
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# Create the main executable
add_executable(agoge_run main.cpp)
target_link_libraries(agoge_run PRIVATE agoge_lib)

# Optionally set RPATH or other linking details here if needed.
