// Example: problems/Sedov.cpp
#include "Sedov.hpp"

#include <cmath>
#include <iostream>

namespace agoge {
namespace problems {

void Sedov::registerParameters(ParameterSystem &params) const {
    std::string prefix = name() + ".";

    params.addDefault(prefix + "radius", "1.0");
    params.addDefault(prefix + "energy", "1.0");
}

void Sedov::initialize(Field3D &Q, const ParameterSystem &params) {
    std::string prefix = name() + ".";

    double radius = params.getDouble(prefix + "radius");
    double energy = params.getDouble(prefix + "energy");
    double edens = energy / (4. / 3. * M_PI * radius * radius * radius);
    // Compute global domain center from Q.global_bbox
    double xMid = (Q.global_bbox.xmin + Q.global_bbox.xmax) / 2.0;
    double yMid = (Q.global_bbox.ymin + Q.global_bbox.ymax) / 2.0;
    double zMid = (Q.global_bbox.zmin + Q.global_bbox.zmax) / 2.0;

    for (int k = 0; k < Q.Nz; ++k) {
        for (int j = 0; j < Q.Ny; ++j) {
            for (int i = 0; i < Q.Nx; ++i) {
                int idx = Q.interiorIndex(i, j, k);

                // Compute cell center coordinates (relative to global domain center)
                double xC = Q.xCenter(i) - xMid;
                double yC = Q.yCenter(j) - yMid;
                double zC = Q.zCenter(k) - zMid;
                double r2 = xC * xC + yC * yC + zC * zC;
                double R2 = radius * radius;

                // Determine if inside the bomb radius
                if (r2 < R2) {
                    Q.rho[idx] = 1.0;
                    Q.E[idx] = edens;
                } else {
                    Q.rho[idx] = 1.0;
                    Q.E[idx] = 0.1;
                }
                // No gravity
                Q.phi[idx] = 0.0;
            }
        }
    }

}

}  // namespace problems
}  // namespace agoge