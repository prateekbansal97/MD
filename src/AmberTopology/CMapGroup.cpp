//
// Created by Prateek Bansal on 12/10/25.
//

#include "../../include/AmberTopology/CMapGroup.h"
#include <vector>
#include <cmath>

double CMapGroup::return_energy_linear(double phi, double psi, const std::vector<double>& grid, int resolution)
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

std::pair<double, double> CMapGroup::return_gradient_linear(double phi, double psi, const std::vector<double>& grid, int resolution)
{
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

    int i_next = (i + 1) % resolution;
    int i_prev = (i - 1 + resolution) % resolution;

    int j_next = (j + 1) % resolution;
    int j_prev = (j - 1 + resolution) % resolution;

    // 4. Look up Energies
    // Grid structure: index = row * width + col
    // We vary i (phi) for the first gradient, keeping j (psi) constant
    double e_phi_next = grid[i_next * resolution + j];
    double e_phi_prev = grid[i_prev * resolution + j];

    // We vary j (psi) for the second gradient, keeping i (phi) constant
    double e_psi_next = grid[i * resolution + j_next];
    double e_psi_prev = grid[i * resolution + j_prev];

    // 5. Calculate Gradients (Central Difference)
    // Slope = Rise / Run = (E_next - E_prev) / (2 * step_size)
    double dE_dphi = (e_phi_next - e_phi_prev) / (2.0 * grid_spacing);
    double dE_dpsi = (e_psi_next - e_psi_prev) / (2.0 * grid_spacing);

    return {dE_dphi, dE_dpsi};
}

double CMapGroup::return_energy_bicubic(double phi, double psi, int resolution, std::vector<double>& coeffs, const std::vector<double>& grid_energies)
{
    // 1. STANDARDIZE ANGLES TO [0, 2*PI]
    double phi_shifted = phi;
    double psi_shifted = psi;

    while (phi_shifted <= -M_PI) phi_shifted += 2.0 * M_PI;
    while (phi_shifted > M_PI)   phi_shifted -= 2.0 * M_PI;
    while (psi_shifted <= -M_PI) psi_shifted += 2.0 * M_PI;
    while (psi_shifted > M_PI)   psi_shifted -= 2.0 * M_PI;

    phi_shifted += M_PI;
    psi_shifted += M_PI;


    double grid_spacing = (2.0 * M_PI) / resolution;


    double grid_pos_phi = phi_shifted / grid_spacing;
    double grid_pos_psi = psi_shifted / grid_spacing;


    int row_grid = static_cast<int>(grid_pos_phi);
    int col_grid = static_cast<int>(grid_pos_psi);


    double u = grid_pos_phi - static_cast<double>(row_grid);
    double v = grid_pos_psi - static_cast<double>(col_grid);


    row_grid = row_grid % resolution;
    col_grid = col_grid % resolution;
    if (row_grid < 0) row_grid += resolution;
    if (col_grid < 0) col_grid += resolution;


    int cmap_stride = resolution * resolution * 16;
    int cmap_start = (parameter_set - 1) * cmap_stride;
    int cell_index = (row_grid * resolution + col_grid);
    int coeff_start_index = cmap_start + (cell_index * 16);


    double u2 = u * u;
    double u3 = u2 * u;
    double v2 = v * v;
    double v3 = v2 * v;

    const double u_pow[4] = {1.0, u, u2, u3};
    const double v_pow[4] = {1.0, v, v2, v3};

    double final_energy = 0.0;
    int k = 0;

    // Loop matches the order we stored: i (u/row) is outer, j (v/col) is inner
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            final_energy += coeffs[coeff_start_index + k] * u_pow[i] * v_pow[j];
            k++;
        }
    }

    return final_energy;
}

//#include <utility> // For std::pair

std::pair<double, double> CMapGroup::return_gradient_bicubic(double phi, double psi, int resolution, std::vector<double>& coeffs, const std::vector<double>& grid_energies)
{
    // 1. STANDARDIZE ANGLES
    double phi_shifted = phi;
    double psi_shifted = psi;

    while (phi_shifted <= -M_PI) phi_shifted += 2.0 * M_PI;
    while (phi_shifted > M_PI)   phi_shifted -= 2.0 * M_PI;
    while (psi_shifted <= -M_PI) psi_shifted += 2.0 * M_PI;
    while (psi_shifted > M_PI)   psi_shifted -= 2.0 * M_PI;

    phi_shifted += M_PI;
    psi_shifted += M_PI;

    // 2. CALCULATE LOCAL COORDINATES (u, v)
    double grid_spacing = (2.0 * M_PI) / resolution;

    // Get float position in grid units
    double grid_pos_phi = phi_shifted / grid_spacing;
    double grid_pos_psi = psi_shifted / grid_spacing;

    // Get integer cell index (unwrapped)
    int row_grid = static_cast<int>(grid_pos_phi);
    int col_grid = static_cast<int>(grid_pos_psi);

    // Calculate u, v BEFORE wrapping indices to avoid boundary discontinuities
    double u = grid_pos_phi - static_cast<double>(row_grid);
    double v = grid_pos_psi - static_cast<double>(col_grid);

    // 3. WRAP INDICES FOR LOOKUP
    row_grid = row_grid % resolution;
    col_grid = col_grid % resolution;
    if (row_grid < 0) row_grid += resolution;
    if (col_grid < 0) col_grid += resolution;

    // 4. LOCATE COEFFICIENTS
    int cmap_stride = resolution * resolution * 16;
    int cmap_start = (parameter_set - 1) * cmap_stride;
    int cell_index = (row_grid * resolution + col_grid);
    int coeff_start_index = cmap_start + (cell_index * 16);

    // 5. PRE-CALCULATE POWERS AND DERIVATIVES
    // Powers for Energy: u^0, u^1, u^2, u^3
    double u2 = u * u;
    double u3 = u2 * u;
    double v2 = v * v;
    double v3 = v2 * v;
    const double u_pow[4] = {1.0, u, u2, u3};
    const double v_pow[4] = {1.0, v, v2, v3};

    // Derivatives of Powers for Gradient: d(u^i)/du
    // d(1)/du = 0
    // d(u)/du = 1
    // d(u^2)/du = 2u
    // d(u^3)/du = 3u^2
    const double du_pow[4] = {0.0, 1.0, 2.0 * u, 3.0 * u2};
    const double dv_pow[4] = {0.0, 1.0, 2.0 * v, 3.0 * v2};

    double dE_du = 0.0;
    double dE_dv = 0.0;
    int k = 0;

    // 6. EVALUATE GRADIENTS
    // dE/du = sum( a_ij * d(u^i)/du * v^j )
    // dE/dv = sum( a_ij * u^i * d(v^j)/dv )

    // Loop order must match storage: i (row/u) outer, j (col/v) inner
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            double coeff = coeffs[coeff_start_index + k];

            // Partial wrt u: derive u term, keep v term constant
            dE_du += coeff * du_pow[i] * v_pow[j];

            // Partial wrt v: keep u term constant, derive v term
            dE_dv += coeff * u_pow[i] * dv_pow[j];

            k++;
        }
    }

    // 7. APPLY CHAIN RULE
    // dE/dPhi = (dE/du) * (du/dPhi)
    // du/dPhi = 1 / grid_spacing
    double inv_spacing = 1.0 / grid_spacing;
    double dE_dphi = dE_du * inv_spacing;
    double dE_dpsi = dE_dv * inv_spacing;

    return {dE_dphi, dE_dpsi};
}












