//
// Created by Prateek Bansal on 12/12/25.
//

#include "../../include/AmberTopology/LennardJones.h"

double LennardJones::CalculateEnergy(double r2, double Aij, double Bij)
{
    const double inv_r2  = 1.0 / r2;
    const double inv_r6  = inv_r2 * inv_r2 * inv_r2;
    const double inv_r12 = inv_r6 * inv_r6;
    return Aij * inv_r12 - Bij * inv_r6;
}