#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>

#include "Field3d.hpp"
#include "Config.hpp"
#include "EulerSolver.hpp" // Include to access transpose functions

// Define IndexOrdering enum since it's not in Field3D.hpp
namespace agoge {
    /// @brief Specifies how 3D data is stored in memory
    enum class IndexOrdering {
        XYZ, // i + j*nx + k*nx*ny (standard row-major, i fastest)
        YXZ, // j + i*ny + k*ny*nx (j fastest)
        ZYX, // k + j*nz + i*nz*ny (k fastest)
        XZY, // i + k*nx + j*nx*nz (i fastest, k second)
        ZXY, // k + i*nz + j*nz*nx (k fastest, i second)
        YZX  // j + k*ny + i*ny*nz (j fastest, k second)
    };
}

/// @brief Helper function to initialize field with position-encoding values
/// Each cell will contain a value of 100*i + 10*j + k
void initializeTestField(agoge::Field3D& field) {
    for (int k = 0; k < field.NzGhost; k++) {
        for (int j = 0; j < field.NyGhost; j++) {
            for (int i = 0; i < field.NxGhost; i++) {
                int idx = field.index(i, j, k);
                double value = 100.0 * i + 10.0 * j + 1.0 * k;
                field.rho[idx] = value;
                field.rhou[idx] = value;
                field.rhov[idx] = value;
                field.rhow[idx] = value;
                field.E[idx] = value;
                field.phi[idx] = value;
            }
        }
    }
}

/// @brief Implementation of getIndex for different orderings
/// This is needed since Field3D doesn't have this method
template<agoge::IndexOrdering Order>
inline int getIndex(const agoge::Field3D& field, int i, int j, int k) {
    const int nx = field.NxGhost;
    const int ny = field.NyGhost;
    const int nz = field.NzGhost;
    
    switch(Order) {
        case agoge::IndexOrdering::XYZ:
            return i + j * nx + k * nx * ny;
        case agoge::IndexOrdering::YXZ:
            return j + i * ny + k * ny * nx;
        case agoge::IndexOrdering::ZYX:
            return k + j * nz + i * nz * ny;
        case agoge::IndexOrdering::XZY:
            return i + k * nx + j * nx * nz;
        case agoge::IndexOrdering::ZXY:
            return k + i * nz + j * nz * nx;
        case agoge::IndexOrdering::YZX:
            return j + k * ny + i * ny * nz;
        default:
            return i + j * nx + k * nx * ny; // Default to XYZ
    }
}

/// @brief Verify that a field contains the expected value at a specific position based on ordering
/// @return true if the field has the expected value, false otherwise
template<agoge::IndexOrdering Order>
bool verifyFieldValue(const agoge::Field3D& field, int i, int j, int k, double expected) {
    // Get the index based on the specified ordering
    int idx = getIndex<Order>(field, i, j, k);
    bool match = std::abs(field.rho[idx] - expected) < 1e-10;
    if (!match) {
        std::cerr << "Value mismatch at (" << i << "," << j << "," << k << "): "
                  << "Expected " << expected << ", got " << field.rho[idx] << std::endl;
    }
    return match;
}

/// @brief Verify multiple field values for a given ordering
template<agoge::IndexOrdering Order>
bool verifyMultipleValues(const agoge::Field3D& field) {
    // Sample points from different regions of the field (interior and ghost)
    std::vector<std::tuple<int, int, int, double>> testPoints = {
        {1, 1, 1, 111.0},                     // Interior point
        {0, 0, 0, 0.0},                       // Corner ghost
        {field.NxGhost-1, 
         field.NyGhost-1, 
         field.NzGhost-1,           // Opposite corner
         100.0*(field.NxGhost-1) + 
         10.0*(field.NyGhost-1) + 
         1.0*(field.NzGhost-1)},
        {agoge::config::Ng, agoge::config::Ng, agoge::config::Ng, // First interior cell
         100.0*agoge::config::Ng + 10.0*agoge::config::Ng + 1.0*agoge::config::Ng}
    };
    
    for (auto& [i, j, k, expected] : testPoints) {
        if (!verifyFieldValue<Order>(field, i, j, k, expected)) {
            return false;
        }
    }
    return true;
}

/// @brief Adapter function to use transpose12 on Field3D data
void transpose12Field(agoge::Field3D& field) {
    const int n1 = field.NxGhost;
    const int n2 = field.NyGhost;
    const int n3 = field.NzGhost;
    
    std::vector<double>* components[] = {&field.rho, &field.rhou, &field.rhov, &field.rhow, &field.E, &field.phi};
    
    for (auto* component : components) {
        std::vector<double> temp(n1 * n2);  // Temporary buffer for one plane
        
        for (int i3 = 0; i3 < n3; ++i3) {
            // Copy current plane to temp buffer with transposed indices
            for (int i2 = 0; i2 < n2; ++i2) {
                for (int i1 = 0; i1 < n1; ++i1) {
                    int srcIdx = i1 + i2 * n1 + i3 * (n1 * n2);
                    int dstIdx = i2 + i1 * n2;
                    temp[dstIdx] = (*component)[srcIdx];
                }
            }
            
            // Copy back to source array with new indexing
            for (int i2 = 0; i2 < n2; ++i2) {
                for (int i1 = 0; i1 < n1; ++i1) {
                    int dstIdx = i1 + i2 * n1 + i3 * (n1 * n2);
                    int tempIdx = i2 + i1 * n2;
                    (*component)[dstIdx] = temp[tempIdx];
                }
            }
        }
    }
}

/// @brief Adapter function to use transpose13 on Field3D data
void transpose13Field(agoge::Field3D& field) {
    const int n1 = field.NxGhost;
    const int n2 = field.NyGhost;
    const int n3 = field.NzGhost;
    
    std::vector<double>* components[] = {&field.rho, &field.rhou, &field.rhov, &field.rhow, &field.E, &field.phi};
    
    for (auto* component : components) {
        std::vector<double> temp(n1 * n3);  // Temporary buffer for one plane
        
        for (int i2 = 0; i2 < n2; ++i2) {
            // Copy current plane to temp buffer with transposed indices
            for (int i3 = 0; i3 < n3; ++i3) {
                for (int i1 = 0; i1 < n1; ++i1) {
                    int srcIdx = i1 + i2 * n1 + i3 * (n1 * n2);
                    int dstIdx = i3 + i1 * n3;
                    temp[dstIdx] = (*component)[srcIdx];
                }
            }
            
            // Copy back to source array with new indexing
            for (int i3 = 0; i3 < n3; ++i3) {
                for (int i1 = 0; i1 < n1; ++i1) {
                    int dstIdx = i3 + i2 * n3 + i1 * (n2 * n3);
                    int tempIdx = i3 + i1 * n3;
                    (*component)[dstIdx] = temp[tempIdx];
                }
            }
        }
    }
}

/// @brief Custom reordering function that uses transpose12 and transpose13
template<agoge::IndexOrdering FromOrder, agoge::IndexOrdering ToOrder>
void reorderField(agoge::Field3D& field) {
    // Determine which transpositions are needed based on FromOrder and ToOrder
    if (FromOrder == agoge::IndexOrdering::XYZ) {
        if (ToOrder == agoge::IndexOrdering::YXZ) {
            transpose12Field(field); // XYZ -> YXZ
        } else if (ToOrder == agoge::IndexOrdering::ZYX) {
            transpose13Field(field); // XYZ -> ZXY
            // For ZXY -> ZYX, swap dims after the first transpose
            agoge::Field3D tempField(field.Nx, field.Ny, field.Nz, field.getBoundingBox(), field.nghost);
            tempField = field;
            transpose12Field(tempField); // ZXY -> ZYX with adjusted dimensions
            field = tempField;
        } else if (ToOrder == agoge::IndexOrdering::XZY) {
            // For XZY, we need a custom approach
            agoge::Field3D tempField(field.Nx, field.Ny, field.Nz, field.getBoundingBox(), field.nghost);
            tempField = field;
            
            // Implement XYZ -> XZY transformation
            const int nx = field.NxGhost;
            const int ny = field.NyGhost;
            const int nz = field.NzGhost;
            std::vector<double>* components[] = {&field.rho, &field.rhou, &field.rhov, &field.rhow, &field.E, &field.phi};
            std::vector<double>* tempComponents[] = {&tempField.rho, &tempField.rhou, &tempField.rhov, &tempField.rhow, &tempField.E, &tempField.phi};
            
            for (int c = 0; c < 6; ++c) {
                for (int i = 0; i < nx; ++i) {
                    for (int j = 0; j < ny; ++j) {
                        for (int k = 0; k < nz; ++k) {
                            int xyzIdx = i + j * nx + k * nx * ny;
                            int xzyIdx = i + k * nx + j * nx * nz;
                            (*components[c])[xzyIdx] = (*tempComponents[c])[xyzIdx];
                        }
                    }
                }
            }
        } else if (ToOrder == agoge::IndexOrdering::ZXY) {
            transpose13Field(field); // XYZ -> ZXY
        } else if (ToOrder == agoge::IndexOrdering::YZX) {
            transpose12Field(field); // XYZ -> YXZ
            transpose13Field(field); // YXZ -> YZX
        }
    } else if (FromOrder == agoge::IndexOrdering::YXZ) {
        if (ToOrder == agoge::IndexOrdering::XYZ) {
            transpose12Field(field); // YXZ -> XYZ
        } else if (ToOrder == agoge::IndexOrdering::ZXY) {
            transpose13Field(field); // YXZ -> ZXY
        }
    } else if (FromOrder == agoge::IndexOrdering::ZXY) {
        if (ToOrder == agoge::IndexOrdering::XYZ) {
            transpose13Field(field); // ZXY -> XYZ
        } else if (ToOrder == agoge::IndexOrdering::YXZ) {
            transpose13Field(field); // ZXY -> XYZ
            transpose12Field(field); // XYZ -> YXZ
        }
    }
}

/// @brief Test direct transformation between two orderings
template<agoge::IndexOrdering FromOrder, agoge::IndexOrdering ToOrder>
bool testDirectTransform() {
    std::cout << "Testing " << static_cast<int>(FromOrder) 
              << " -> " << static_cast<int>(ToOrder) << " transform..." << std::endl;
    
    // Create a small test field with appropriate dimensions
    const int nx = 32, ny = 32, nz = 32;
    const int nghost = 2;
    agoge::BoundingBox box = {0.0, 1.0, 0.0, 1.0, 0.0, 1.0};
    agoge::Field3D field(nx, ny, nz, box, nghost);
    
    // Initialize with position-encoding values
    initializeTestField(field);
    
    // Verify initial values
    if (!verifyMultipleValues<FromOrder>(field)) {
        std::cerr << "Initial field verification failed" << std::endl;
        return false;
    }
    
    // Perform the reordering using our new function
    reorderField<FromOrder, ToOrder>(field);
    
    // Verify reordered values
    bool success = verifyMultipleValues<ToOrder>(field);
    
    if (success) {
        std::cout << "Transform " << static_cast<int>(FromOrder) 
                  << " -> " << static_cast<int>(ToOrder) << " successful!" << std::endl;
    } else {
        std::cerr << "Transform " << static_cast<int>(FromOrder) 
                  << " -> " << static_cast<int>(ToOrder) << " failed!" << std::endl;
    }
    
    return success;
}

/// @brief Test round-trip transformation (FromOrder -> ToOrder -> FromOrder)
template<agoge::IndexOrdering FromOrder, agoge::IndexOrdering ToOrder>
bool testRoundTripTransform() {
    std::cout << "Testing round-trip " << static_cast<int>(FromOrder) 
              << " -> " << static_cast<int>(ToOrder) 
              << " -> " << static_cast<int>(FromOrder) << "..." << std::endl;
    
    // Create a small test field with appropriate dimensions
    const int nx = 32, ny = 32, nz = 32;
    const int nghost = 2;
    agoge::BoundingBox box = {0.0, 1.0, 0.0, 1.0, 0.0, 1.0};
    agoge::Field3D field(nx, ny, nz, box, nghost);
    
    // Save original data for later comparison
    initializeTestField(field);
    std::vector<double> originalData = field.rho; // Save a copy
    
    // First transformation
    reorderField<FromOrder, ToOrder>(field);
    
    // Second transformation (back to original)
    reorderField<ToOrder, FromOrder>(field);
    
    // Verify data is restored to original
    bool success = true;
    for (size_t i = 0; i < field.rho.size(); i++) {
        if (std::abs(field.rho[i] - originalData[i]) > 1e-10) {
            std::cerr << "Round-trip mismatch at index " << i 
                      << ": Expected " << originalData[i] 
                      << ", got " << field.rho[i] << std::endl;
            success = false;
            break;
        }
    }
    
    if (success) {
        std::cout << "Round-trip " << static_cast<int>(FromOrder) 
                  << " -> " << static_cast<int>(ToOrder) 
                  << " -> " << static_cast<int>(FromOrder) << " successful!" << std::endl;
    } else {
        std::cerr << "Round-trip " << static_cast<int>(FromOrder) 
                  << " -> " << static_cast<int>(ToOrder) 
                  << " -> " << static_cast<int>(FromOrder) << " failed!" << std::endl;
    }
    
    return success;
}

int main() {
    bool allTestsPassed = true;
    
    // Test direct transformations between all ordering pairs
    allTestsPassed &= testDirectTransform<agoge::IndexOrdering::XYZ, agoge::IndexOrdering::YXZ>();
    allTestsPassed &= testDirectTransform<agoge::IndexOrdering::XYZ, agoge::IndexOrdering::ZYX>();
    allTestsPassed &= testDirectTransform<agoge::IndexOrdering::XYZ, agoge::IndexOrdering::XZY>();
    allTestsPassed &= testDirectTransform<agoge::IndexOrdering::XYZ, agoge::IndexOrdering::ZXY>();
    allTestsPassed &= testDirectTransform<agoge::IndexOrdering::XYZ, agoge::IndexOrdering::YZX>();
    
    allTestsPassed &= testDirectTransform<agoge::IndexOrdering::YXZ, agoge::IndexOrdering::XYZ>();
    allTestsPassed &= testDirectTransform<agoge::IndexOrdering::YXZ, agoge::IndexOrdering::ZXY>();
    
    allTestsPassed &= testDirectTransform<agoge::IndexOrdering::ZXY, agoge::IndexOrdering::XYZ>();
    allTestsPassed &= testDirectTransform<agoge::IndexOrdering::ZXY, agoge::IndexOrdering::YXZ>();
    
    // Test round-trip transformations for key ordering combinations
    allTestsPassed &= testRoundTripTransform<agoge::IndexOrdering::XYZ, agoge::IndexOrdering::YXZ>();
    allTestsPassed &= testRoundTripTransform<agoge::IndexOrdering::XYZ, agoge::IndexOrdering::ZXY>();
    allTestsPassed &= testRoundTripTransform<agoge::IndexOrdering::YXZ, agoge::IndexOrdering::ZXY>();
    
    // More complex scenario: test XYZ -> YXZ -> ZXY -> XYZ
    std::cout << "Testing complex transformation sequence XYZ -> YXZ -> ZXY -> XYZ..." << std::endl;
    {
        const int nx = 32, ny = 32, nz = 32;
        const int nghost = 2;
        agoge::BoundingBox box = {0.0, 1.0, 0.0, 1.0, 0.0, 1.0};
        agoge::Field3D field(nx, ny, nz, box, nghost);
        
        initializeTestField(field);
        std::vector<double> originalData = field.rho;
        
        reorderField<agoge::IndexOrdering::XYZ, agoge::IndexOrdering::YXZ>(field);
        reorderField<agoge::IndexOrdering::YXZ, agoge::IndexOrdering::ZXY>(field);
        reorderField<agoge::IndexOrdering::ZXY, agoge::IndexOrdering::XYZ>(field);
        
        bool success = true;
        for (size_t i = 0; i < field.rho.size(); i++) {
            if (std::abs(field.rho[i] - originalData[i]) > 1e-10) {
                std::cerr << "Complex transformation mismatch at index " << i << std::endl;
                success = false;
                break;
            }
        }
        
        if (success) {
            std::cout << "Complex transformation successful!" << std::endl;
        } else {
            std::cerr << "Complex transformation failed!" << std::endl;
            allTestsPassed = false;
        }
    }
    
    if (allTestsPassed) {
        std::cout << "All tests passed!" << std::endl;
        return 0;
    } else {
        std::cerr << "Some tests failed!" << std::endl;
        return 1;
    }
}
