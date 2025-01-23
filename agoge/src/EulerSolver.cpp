#include "EulerSolver.hpp"
#include "Field3d.hpp"
#include "Config.hpp"

// NEW:
#include "PerformanceMonitor.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <array>

namespace agoge
{
    namespace euler
    {

        static inline double pressure(double rho, double rhou, double rhov,
                                      double rhow, double E)
        {
            double gamma_gas = agoge::config::gamma_gas;
            double u = rhou / rho;
            double v = rhov / rho;
            double w = rhow / rho;
            double kinetic = 0.5 * rho * (u * u + v * v + w * w);
            double p = (gamma_gas - 1.0) * (E - kinetic);
            return p;
        }

        static inline void computeGravityAccel(const Field3D &gravField,
                                               int i, int j, int k,
                                               double &gx, double &gy, double &gz)
        {
            int Nx = gravField.Nx;
            int Ny = gravField.Ny;
            int Nz = gravField.Nz;

            double dx = gravField.dx;
            double dy = gravField.dy;
            double dz = gravField.dz;

            int iL = (i == 0) ? (Nx - 1) : (i - 1);
            int iR = (i == Nx - 1) ? 0 : (i + 1);
            int jL = (j == 0) ? (Ny - 1) : (j - 1);
            int jR = (j == Ny - 1) ? 0 : (j + 1);
            int kL = (k == 0) ? (Nz - 1) : (k - 1);
            int kR = (k == Nz - 1) ? 0 : (k + 1);

            int idxC = gravField.index(i, j, k);
            int idxXm = gravField.index(iL, j, k);
            int idxXp = gravField.index(iR, j, k);
            int idxYm = gravField.index(i, jL, k);
            int idxYp = gravField.index(i, jR, k);
            int idxZm = gravField.index(i, j, kL);
            int idxZp = gravField.index(i, j, kR);

            double phiXm = gravField.phi[idxXm];
            double phiXp = gravField.phi[idxXp];
            double phiYm = gravField.phi[idxYm];
            double phiYp = gravField.phi[idxYp];
            double phiZm = gravField.phi[idxZm];
            double phiZp = gravField.phi[idxZp];

            double dphidx = (phiXp - phiXm) / (2.0 * dx);
            double dphidy = (phiYp - phiYm) / (2.0 * dy);
            double dphidz = (phiZp - phiZm) / (2.0 * dz);

            gx = -dphidx;
            gy = -dphidy;
            gz = -dphidz;
        }

        void computeL(const Field3D &Q, Field3D &LQ, const Field3D *gravField)
        {
            // START TIMER for computeL
            agoge::PerformanceMonitor::instance().startTimer("computeL");

            std::fill(LQ.rho.begin(), LQ.rho.end(), 0.0);
            std::fill(LQ.rhou.begin(), LQ.rhou.end(), 0.0);
            std::fill(LQ.rhov.begin(), LQ.rhov.end(), 0.0);
            std::fill(LQ.rhow.begin(), LQ.rhow.end(), 0.0);
            std::fill(LQ.E.begin(), LQ.E.end(), 0.0);
            std::fill(LQ.phi.begin(), LQ.phi.end(), 0.0);

            int Nx = Q.Nx;
            int Ny = Q.Ny;
            int Nz = Q.Nz;

            // Define flux functions with exactly 5 elements
            auto fluxX = [&](double r, double ru, double rv, double rw, double e) -> std::array<double, 5>
            {
                double pcell = pressure(r, ru, rv, rw, e);
                double u = ru / r;
                return std::array<double, 5>{
                    ru,             // Mass flux
                    ru * u + pcell, // Momentum flux x
                    ru * rv / r,    // Momentum flux y
                    ru * rw / r,    // Momentum flux z
                    (e + pcell) * u // Energy flux
                };
            };

            auto fluxY = [&](double r, double ru, double rv, double rw, double e) -> std::array<double, 5>
            {
                double pcell = pressure(r, ru, rv, rw, e);
                double v = rv / r;
                return std::array<double, 5>{
                    rv,             // Mass flux
                    ru * v,         // Momentum flux x
                    rv * v + pcell, // Momentum flux y
                    rv * rw / r,    // Momentum flux z
                    (e + pcell) * v // Energy flux
                };
            };

            auto fluxZ = [&](double r, double ru, double rv, double rw, double e) -> std::array<double, 5>
            {
                double pcell = pressure(r, ru, rv, rw, e);
                double w = rw / r;
                return std::array<double, 5>{
                    rw,             // Mass flux
                    ru * w,         // Momentum flux x
                    rv * w,         // Momentum flux y
                    rw * w + pcell, // Momentum flux z
                    (e + pcell) * w // Energy flux
                };
            };

            // ... Rest of your computeL implementation ...

            // STOP TIMER for computeL
            agoge::PerformanceMonitor::instance().stopTimer("computeL");
        }

        void runRK2(Field3D &Q, double dt)
        {
            agoge::PerformanceMonitor::instance().startTimer("runRK2");

            Field3D Qtemp(Q.Nx, Q.Ny, Q.Nz, Q.dx, Q.dy, Q.dz);
            Field3D LQ(Q.Nx, Q.Ny, Q.Nz, Q.dx, Q.dy, Q.dz);
            Field3D LQtemp(Q.Nx, Q.Ny, Q.Nz, Q.dx, Q.dy, Q.dz);

            // Stage 1
            computeL(Q, LQ, nullptr); // If you store gravity in Q, pass &Q
            for (size_t n = 0; n < Q.rho.size(); ++n)
            {
                Qtemp.rho[n] = Q.rho[n] + dt * LQ.rho[n];
                Qtemp.rhou[n] = Q.rhou[n] + dt * LQ.rhou[n];
                Qtemp.rhov[n] = Q.rhov[n] + dt * LQ.rhov[n];
                Qtemp.rhow[n] = Q.rhow[n] + dt * LQ.rhow[n];
                Qtemp.E[n] = Q.E[n] + dt * LQ.E[n];
            }

            // Stage 2
            computeL(Qtemp, LQtemp, nullptr); // If you store gravity in Q, pass &Qtemp
            for (size_t n = 0; n < Q.rho.size(); ++n)
            {
                Q.rho[n] = 0.5 * (Q.rho[n] + Qtemp.rho[n] + dt * LQtemp.rho[n]);
                Q.rhou[n] = 0.5 * (Q.rhou[n] + Qtemp.rhou[n] + dt * LQtemp.rhou[n]);
                Q.rhov[n] = 0.5 * (Q.rhov[n] + Qtemp.rhov[n] + dt * LQtemp.rhov[n]);
                Q.rhow[n] = 0.5 * (Q.rhow[n] + Qtemp.rhow[n] + dt * LQtemp.rhow[n]);
                Q.E[n] = 0.5 * (Q.E[n] + Qtemp.E[n] + dt * LQtemp.E[n]);
            }

            agoge::PerformanceMonitor::instance().stopTimer("runRK2");
        }

        double computeTimeStep(const Field3D &Q, double cfl)
        {
            agoge::PerformanceMonitor::instance().startTimer("computeTimeStep");
            double maxSpeed = 0.0;
            double gamma_gas = agoge::config::gamma_gas;

            int Nx = Q.Nx;
            int Ny = Q.Ny;
            int Nz = Q.Nz;

            for (int k = 0; k < Nz; ++k)
            {
                for (int j = 0; j < Ny; ++j)
                {
                    for (int i = 0; i < Nx; ++i)
                    {
                        int idx = Q.index(i, j, k);
                        double r = Q.rho[idx];
                        double ru = Q.rhou[idx];
                        double rv = Q.rhov[idx];
                        double rw = Q.rhow[idx];
                        double e = Q.E[idx];

                        if (r <= 0.0)
                            continue;

                        double u = ru / r;
                        double v = rv / r;
                        double w = rw / r;
                        double speed = std::sqrt(u * u + v * v + w * w);

                        double p = pressure(r, ru, rv, rw, e);
                        if (p < 0.0)
                            continue;

                        double a = std::sqrt(gamma_gas * p / r);
                        double localWave = speed + a;
                        if (localWave > maxSpeed)
                        {
                            maxSpeed = localWave;
                        }
                    }
                }
            }

            double minDx = std::min({Q.dx, Q.dy, Q.dz});
            double dt = 1.e20;
            if (maxSpeed > 1.e-14)
            {
                dt = cfl * (minDx / maxSpeed);
            }

            agoge::PerformanceMonitor::instance().stopTimer("computeTimeStep");
            return dt;
        }

    } // namespace euler
} // namespace agoge
