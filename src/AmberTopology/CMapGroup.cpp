//
// Created by Prateek Bansal on 12/10/25.
//

#include "../../include/AmberTopology/CMapGroup.h"
#include <vector>
#include <cmath>

double CMapGroup::return_energy(double phi, double psi, const std::vector<double>& grid, int resolution)
{
    // 1. Shift angles from [-PI, PI] to [0, 2*PI]
    double phi_shifted = phi;
    double psi_shifted = psi;

    // Standardize range to [-PI, PI] first, then shift
    while (phi_shifted <= -M_PI) phi_shifted += 2.0 * M_PI;
    while (phi_shifted > M_PI)   phi_shifted -= 2.0 * M_PI;
    while (psi_shifted <= -M_PI) psi_shifted += 2.0 * M_PI;
    while (psi_shifted > M_PI)   psi_shifted -= 2.0 * M_PI;

    // Shift to [0, 2*PI] for indexing
    phi_shifted += M_PI;
    psi_shifted += M_PI;

    double grid_spacing = (2.0 * M_PI) / resolution;
    int i = static_cast<int>(phi_shifted / grid_spacing) % resolution;
    int j = static_cast<int>(psi_shifted / grid_spacing) % resolution;

    if (i < 0) i += resolution;
    if (j < 0) j += resolution;

    return grid[i * resolution + j];
}
