//
// Created by Prateek Bansal on 12/9/25.
//
#include "../../include/AmberTopology/topology.h"
#include "../../include/System/System.h"
#include <vector>
#include <cmath>
#include <algorithm>

void System::init() {
    const std::vector<double>& coordinates = topology.get_coordinates();

    double bond_energies = 0;
    for (auto& bond: topology.get_harmonic_bonds())
    {
        int atomAIndex = bond.get_atomA_index();
        int atomBIndex = bond.get_atomB_index();
        const double x1 = coordinates[3*atomAIndex], y1 = coordinates[3*atomAIndex + 1], z1 = coordinates[3*atomAIndex + 2];
        const double x2 = coordinates[3*atomBIndex], y2 = coordinates[3*atomBIndex + 1], z2 = coordinates[3*atomBIndex + 2];
        double dist = distance(x1, y1, z1, x2, y2, z2);
        bond.set_distance(dist);

        double energy = bond.return_energy(dist);
        bond.set_energy(energy);
        bond_energies += energy;
    }
    double angle_energy = 0;
    for (auto& ang: topology.get_harmonic_angles())
    {
        int atomAIndex = ang.get_atomA_index();
        int atomBIndex = ang.get_atomB_index();
        int atomCIndex = ang.get_atomC_index();
        const double x1 = coordinates[3*atomAIndex], y1 = coordinates[3*atomAIndex + 1], z1 = coordinates[3*atomAIndex + 2];
        const double x2 = coordinates[3*atomBIndex], y2 = coordinates[3*atomBIndex + 1], z2 = coordinates[3*atomBIndex + 2];
        const double x3 = coordinates[3*atomCIndex], y3 = coordinates[3*atomCIndex + 1], z3 = coordinates[3*atomCIndex + 2];
        double angle_m = angle(x1, y1, z1, x2, y2, z2, x3, y3, z3);
        ang.set_angle(angle_m);

        double energy = ang.return_energy(angle_m);
        ang.set_energy(energy);
        angle_energy += energy;
    }

    double dihedral_energy = 0;
    for (auto& dih: topology.get_cosine_dihedrals())
    {
        int atomAIndex = dih.get_atomA_index();
        int atomBIndex = dih.get_atomB_index();
        int atomCIndex = dih.get_atomC_index();
        int atomDIndex = dih.get_atomD_index();
        const double x1 = coordinates[3*atomAIndex], y1 = coordinates[3*atomAIndex + 1], z1 = coordinates[3*atomAIndex + 2];
        const double x2 = coordinates[3*atomBIndex], y2 = coordinates[3*atomBIndex + 1], z2 = coordinates[3*atomBIndex + 2];
        const double x3 = coordinates[3*atomCIndex], y3 = coordinates[3*atomCIndex + 1], z3 = coordinates[3*atomCIndex + 2];
        const double x4 = coordinates[3*atomDIndex], y4 = coordinates[3*atomDIndex + 1], z4 = coordinates[3*atomDIndex + 2];
        double angle_m = dihedral(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4);
        dih.set_cosine_dihedral(angle_m);
        double energy = dih.return_energy(angle_m);
        dih.set_energy(energy);
        dihedral_energy += energy;
    }

    double urey_bradley_energy = 0;
    for (auto& ub : topology.get_harmonic_UBs())
    {
        int a1 = ub.get_atomA_index();
        int a2 = ub.get_atomB_index();

        double dist = distance(
                coordinates[3*a1], coordinates[3*a1+1], coordinates[3*a1+2],
                coordinates[3*a2], coordinates[3*a2+1], coordinates[3*a2+2]
        );

         ub.set_distance_value(dist);
         double energy = ub.return_energy(dist);
         ub.set_energy(energy);
         urey_bradley_energy += energy;
//         std::cout << "UB Distance: " << dist << "\n";
    }

    double improper_energy = 0;
    for (auto& imp : topology.get_harmonic_impropers())
    {
        int a1 = imp.get_atomA_index();
        int a2 = imp.get_atomB_index();
        int a3 = imp.get_atomC_index();
        int a4 = imp.get_atomD_index();

        double angle_m = dihedral(
                coordinates[3*a1], coordinates[3*a1+1], coordinates[3*a1+2],
                coordinates[3*a2], coordinates[3*a2+1], coordinates[3*a2+2],
                coordinates[3*a3], coordinates[3*a3+1], coordinates[3*a3+2],
                coordinates[3*a4], coordinates[3*a4+1], coordinates[3*a4+2]
        );

        // Assuming you added a 'current_angle' member or setter to HarmonicImproper
         imp.set_imp_dihedral(angle_m);
         double energy = imp.return_energy(angle_m);
         imp.set_energy(energy);
         improper_energy += energy;

        // Debug print (Optional)
//         std::cout << "Improper Angle: " << angle_m*180/M_PI << "\n";
    }

    double cmap_energy = 0.0;
    for (auto& cmap : topology.get_cmaps())
    {
        int a1 = cmap.get_atomA_index();
        int a2 = cmap.get_atomB_index();
        int a3 = cmap.get_atomC_index();
        int a4 = cmap.get_atomD_index();
        int a5 = cmap.get_atomE_index();

        // Calculate First Torsion (A-B-C-D)
        double phi = dihedral(
                coordinates[3*a1], coordinates[3*a1+1], coordinates[3*a1+2],
                coordinates[3*a2], coordinates[3*a2+1], coordinates[3*a2+2],
                coordinates[3*a3], coordinates[3*a3+1], coordinates[3*a3+2],
                coordinates[3*a4], coordinates[3*a4+1], coordinates[3*a4+2]
        );

        // Calculate Second Torsion (B-C-D-E)
        double psi = dihedral(
                coordinates[3*a2], coordinates[3*a2+1], coordinates[3*a2+2],
                coordinates[3*a3], coordinates[3*a3+1], coordinates[3*a3+2],
                coordinates[3*a4], coordinates[3*a4+1], coordinates[3*a4+2],
                coordinates[3*a5], coordinates[3*a5+1], coordinates[3*a5+2]
        );

        int set_index = cmap.get_parameter_set();
        const std::vector<int>& resolutions = topology.get_charmm_cmap_resolutions();
        int resolution = resolutions[set_index - 1];
        const std::vector<double>& full_grid = topology.get_cmap_grid_data(set_index);


        // Pass the OFFSET and the FULL GRID to your return_energy function
        double energy = cmap.return_energy(phi, psi, full_grid, resolution);

        cmap.set_angles(phi, psi);
        cmap.set_energy(energy);
        cmap_energy += energy;
        // Debug
//         std::cout << "CMAP Angles: " << phi*180/M_PI << ", " << psi*180/M_PI << "\n";
    }
    std::cout << "\n Bond:" << bond_energies << " Angle: " << angle_energy << " CosineDihedral: " << dihedral_energy <<
    " Urey-Bradley: " << urey_bradley_energy << " Impropers: " << improper_energy << " CMAP energy: " << cmap_energy << "\n";
}

double System::distance(double x1, double y1, double z1, double x2, double y2, double z2)
{
    const double dx = x1 - x2;
    const double dy = y1 - y2;
    const double dz = z1 - z2;
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

double System::angle(double x1, double y1, double z1, double x2, double y2, double z2,  double x3, double y3, double z3)
{
    const double dx1 = x1 - x2;
    const double dy1 = y1 - y2;
    const double dz1 = z1 - z2;

    const double dx2 = x3 - x2;
    const double dy2 = y3 - y2;
    const double dz2 = z3 - z2;

    const double dot = dx1*dx2 + dy1*dy2 + dz1*dz2;
    const double mag1_sq = dx1*dx1 + dy1*dy1 + dz1*dz1;
    const double mag2_sq = dx2*dx2 + dy2*dy2 + dz2*dz2;
    const double eps = 1e-14;
    if (mag1_sq < eps || mag2_sq < eps) {
        return std::numeric_limits<double>::quiet_NaN(); // or 0.0, or throw
    }
    const double inv_mag = 1.0 / std::sqrt(mag1_sq * mag2_sq);
    double c = dot * inv_mag;
    c = std::clamp(c, -1.0, 1.0);
    return std::acos(c);
}

double System::dihedral(double x1, double y1, double z1, double x2, double y2, double z2,
                double x3, double y3, double z3, double x4, double y4, double z4)
{
    // Calculating P1, P2, P3 as 2 - 1, 3 - 2, 4 - 3
    const double dx1 = x2 - x1;
    const double dy1 = y2 - y1;
    const double dz1 = z2 - z1;
    const double dx2 = x3 - x2;
    const double dy2 = y3 - y2;
    const double dz2 = z3 - z2;

    const double dx3 = x4 - x3;
    const double dy3 = y4 - y3;
    const double dz3 = z4 - z3;

    // Calculating n1 = P1 x P2
    double n1x = dy1*dz2 - dz1*dy2;
    double n1y = dx2*dz1 - dx1*dz2; //-1 * (dx1*dz2 - dx2*dz1);
    double n1z = dx1*dy2 - dy1*dx2;

    // Calculating n2 = P2 x P3
    double n2x = dy2*dz3 - dz2*dy3;
    double n2y = dx3*dz2 - dx2*dz3;//-1 * (dx2*dz3 - dx3*dz2);
    double n2z = dx2*dy3 - dy2*dx3;

    const double magn1 = n1x*n1x + n1y*n1y + n1z*n1z;
    const double magn2 = n2x*n2x + n2y*n2y + n2z*n2z;
    const double magb1 = std::sqrt(dx2*dx2+dy2*dy2+dz2*dz2);

    const double eps = 1e-14;
    const double eps2 = 1e-28;

    if (magn1 < eps2 || magn2 < eps2 || magb1 < eps) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double x = n1x*n2x + n1y*n2y + n1z*n2z; // n1 dot n2

    double val = dx1*n2x + dy1*n2y + dz1*n2z; // n2 dot b1
    const double y = val*magb1;
    return std::atan2(y, x);
}