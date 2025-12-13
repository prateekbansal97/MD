//
// Created by Prateek Bansal on 12/12/25.
//

#include "../../include/AmberTopology/LennardJones.h"
#include <cmath>

double LennardJones::CalculateEnergy(double r2, double Aij, double Bij)
{
    const double inv_r2  = 1.0 / r2;
    const double inv_r6  = inv_r2 * inv_r2 * inv_r2;
    const double inv_r12 = inv_r6 * inv_r6;
    return Aij * inv_r12 - Bij * inv_r6;
}

double LennardJones::CalculateGradient(double r2, double Aij, double Bij)
{
    const double inv_r2  = 1.0 / r2;                // 1/r^2
    const double inv_r4  = inv_r2 * inv_r2;         // 1/r^4
    const double inv_r6  = inv_r4 * inv_r2;         // 1/r^6
    const double inv_r8  = inv_r4 * inv_r4;         // 1/r^8
    const double inv_r12 = inv_r6 * inv_r6;         // 1/r^12
    const double inv_r14 = inv_r12 * inv_r2;        // 1/r^14

    return 12.0 * Aij * inv_r14 - 6.0 * Bij * inv_r8;
}