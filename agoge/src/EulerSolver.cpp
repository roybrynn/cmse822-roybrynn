
//------------------------------------------------------------
// Helper functions
//------------------------------------------------------------

// Pressure
inline double pressure(double rho, double rhou, double rhov, double rhow, double E)
{
    double u = rhou / rho;
    double v = rhov / rho;
    double w = rhow / rho;
    double kinetic = 0.5 * rho * (u*u + v*v + w*w);
    double p = (gamma_gas - 1.0) * (E - kinetic);
    return p;
}

// Maximum wave speed (for CFL)
inline double maxWaveSpeed(double rho, double rhou, double rhov, double rhow, double E)
{
    double p = pressure(rho, rhou, rhov, rhow, E);
    double u = rhou / rho;
    double v = rhov / rho;
    double w = rhow / rho;
    double a = std::sqrt(gamma_gas * p / rho); // speed of sound
    double speed = std::sqrt(u*u + v*v + w*w);
    return speed + a;
}

// x-direction flux
inline void fluxX(double rho, double rhou, double rhov, double rhow, double E,
                  double& Fx0, double& Fx1, double& Fx2, double& Fx3, double& Fx4)
{
    double p = pressure(rho, rhou, rhov, rhow, E);
    double u = rhou / rho;
    Fx0 = rhou;
    Fx1 = rhou*u + p;         // rho*u^2 + p
    Fx2 = rhou*(rhov / rho);  // rho*u*v
    Fx3 = rhou*(rhow / rho);  // rho*u*w
    Fx4 = (E + p)*u;
}

// y-direction flux
inline void fluxY(double rho, double rhou, double rhov, double rhow, double E,
                  double& Fy0, double& Fy1, double& Fy2, double& Fy3, double& Fy4)
{
    double p = pressure(rho, rhou, rhov, rhow, E);
    double v = rhov / rho;
    Fy0 = rhov;
    Fy1 = rhov*(rhou / rho);  // rho*v*u
    Fy2 = rhov*v + p;         // rho*v^2 + p
    Fy3 = rhov*(rhow / rho);  // rho*v*w
    Fy4 = (E + p)*v;
}

// z-direction flux
inline void fluxZ(double rho, double rhou, double rhov, double rhow, double E,
                  double& Fz0, double& Fz1, double& Fz2, double& Fz3, double& Fz4)
{
    double p = pressure(rho, rhou, rhov, rhow, E);
    double w = rhow / rho;
    Fz0 = rhow;
    Fz1 = rhow*(rhou / rho);  // rho*w*u
    Fz2 = rhow*(rhov / rho);  // rho*w*v
    Fz3 = rhow*w + p;         // rho*w^2 + p
    Fz4 = (E + p)*w;
}

//------------------------------------------------------------
// Compute L(Q) = - ( dF/dx + dF/dy + dF/dz ) using central diffs
//------------------------------------------------------------
void computeL(const Field3D& Q,
              std::vector<double>& L_rho,
              std::vector<double>& L_rhou,
              std::vector<double>& L_rhov,
              std::vector<double>& L_rhow,
              std::vector<double>& L_E)
{
    const int Nx_ = Q.Nx_;
    const int Ny_ = Q.Ny_;
    const int Nz_ = Q.Nz_;
    const double dx_ = Q.dx_;
    const double dy_ = Q.dy_;
    const double dz_ = Q.dz_;

    std::fill(L_rho.begin(),  L_rho.end(),  0.0);
    std::fill(L_rhou.begin(), L_rhou.end(), 0.0);
    std::fill(L_rhov.begin(), L_rhov.end(), 0.0);
    std::fill(L_rhow.begin(), L_rhow.end(), 0.0);
    std::fill(L_E.begin(),    L_E.end(),    0.0);

    double FxL[5], FxR[5];
    double FyL[5], FyR[5];
    double FzL[5], FzR[5];

    // #pragma omp parallel for collapse(3)  // Optionally parallelize
    for(int k = 0; k < Nz_; k++) {
        for(int j = 0; j < Ny_; j++) {
            for(int i = 0; i < Nx_; i++) {

                // Periodic in x, y, z
                int iL = (i == 0)      ? (Nx_ - 1) : (i - 1);
                int iR = (i == Nx_ - 1)? 0         : (i + 1);

                int jL = (j == 0)      ? (Ny_ - 1) : (j - 1);
                int jR = (j == Ny_ - 1)? 0         : (j + 1);

                int kL = (k == 0)      ? (Nz_ - 1) : (k - 1);
                int kR = (k == Nz_ - 1)? 0         : (k + 1);

                // Current cell
                int idxC = Q.index(i,j,k);

                // Neighbor indices
                int idxXm = Q.index(iL,j,k), idxXp = Q.index(iR,j,k);
                int idxYm = Q.index(i,jL,k), idxYp = Q.index(i,jR,k);
                int idxZm = Q.index(i,j,kL), idxZp = Q.index(i,j,kR);

                // Load neighbor data
                fluxX(Q.rho [idxXm], Q.rhou[idxXm], Q.rhov[idxXm], Q.rhow[idxXm], Q.E[idxXm],
                      FxL[0], FxL[1], FxL[2], FxL[3], FxL[4]);
                fluxX(Q.rho [idxXp], Q.rhou[idxXp], Q.rhov[idxXp], Q.rhow[idxXp], Q.E[idxXp],
                      FxR[0], FxR[1], FxR[2], FxR[3], FxR[4]);

                fluxY(Q.rho [idxYm], Q.rhou[idxYm], Q.rhov[idxYm], Q.rhow[idxYm], Q.E[idxYm],
                      FyL[0], FyL[1], FyL[2], FyL[3], FyL[4]);
                fluxY(Q.rho [idxYp], Q.rhou[idxYp], Q.rhov[idxYp], Q.rhow[idxYp], Q.E[idxYp],
                      FyR[0], FyR[1], FyR[2], FyR[3], FyR[4]);

                fluxZ(Q.rho [idxZm], Q.rhou[idxZm], Q.rhov[idxZm], Q.rhow[idxZm], Q.E[idxZm],
                      FzL[0], FzL[1], FzL[2], FzL[3], FzL[4]);
                fluxZ(Q.rho [idxZp], Q.rhou[idxZp], Q.rhov[idxZp], Q.rhow[idxZp], Q.E[idxZp],
                      FzR[0], FzR[1], FzR[2], FzR[3], FzR[4]);

                double inv2dx = 1.0/(2.0*dx_);
                double inv2dy = 1.0/(2.0*dy_);
                double inv2dz = 1.0/(2.0*dz_);

                // Differences
                double dFxdx0 = (FxR[0] - FxL[0]) * inv2dx;
                double dFxdx1 = (FxR[1] - FxL[1]) * inv2dx;
                double dFxdx2 = (FxR[2] - FxL[2]) * inv2dx;
                double dFxdx3 = (FxR[3] - FxL[3]) * inv2dx;
                double dFxdx4 = (FxR[4] - FxL[4]) * inv2dx;

                double dFydy0 = (FyR[0] - FyL[0]) * inv2dy;
                double dFydy1 = (FyR[1] - FyL[1]) * inv2dy;
                double dFydy2 = (FyR[2] - FyL[2]) * inv2dy;
                double dFydy3 = (FyR[3] - FyL[3]) * inv2dy;
                double dFydy4 = (FyR[4] - FyL[4]) * inv2dy;

                double dFzdz0 = (FzR[0] - FzL[0]) * inv2dz;
                double dFzdz1 = (FzR[1] - FzL[1]) * inv2dz;
                double dFzdz2 = (FzR[2] - FzL[2]) * inv2dz;
                double dFzdz3 = (FzR[3] - FzL[3]) * inv2dz;
                double dFzdz4 = (FzR[4] - FzL[4]) * inv2dz;

                // Accumulate in L arrays
                L_rho [idxC] = - (dFxdx0 + dFydy0 + dFzdz0);
                L_rhou[idxC] = - (dFxdx1 + dFydy1 + dFzdz1);
                L_rhov[idxC] = - (dFxdx2 + dFydy2 + dFzdz2);
                L_rhow[idxC] = - (dFxdx3 + dFydy3 + dFzdz3);
                L_E   [idxC] = - (dFxdx4 + dFydy4 + dFzdz4);

                // Compute partial derivatives of phi (central difference)
                double dphidx = ( phi[idx(i+1,j,k)] - phi[idx(i-1,j,k)] ) / (2.0*dx );
                double dphidy = ( phi[idx(i,j+1,k)] - phi[idx(i,j-1,k)] ) / (2.0*dy );
                double dphidz = ( phi[idx(i,j,k+1)] - phi[idx(i,j,k-1)] ) / (2.0*dz );

                // Gravitational acceleration components
                double gx = -dphidx;
                double gy = -dphidy;
                double gz = -dphidz;

                // Then add to L_rhou, L_rhov, L_rhow
                L_rhou[idxC] += Q.rho[idxC] * gx;
                L_rhov[idxC] += Q.rho[idxC] * gy;
                L_rhow[idxC] += Q.rho[idxC] * gz;

                // And to L_E:  ( -rho * u dot grad(phi) ) = rho * (u dot g)
                double ux = Q.rhou[idxC]/Q.rho[idxC];
                double uy = Q.rhov[idxC]/Q.rho[idxC];
                double uz = Q.rhow[idxC]/Q.rho[idxC];

                L_E[idxC] += Q.rho[idxC] * (ux*gx + uy*gy + uz*gz);
            }
        }
    }
}
