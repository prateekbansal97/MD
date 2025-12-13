//
// Created by Prateek Bansal on 12/9/25.
//
#include "../../include/AmberTopology/topology.h"
#include "../../include/System/System.h"
#include "../../include/AmberTopology/LennardJones.h"
#include "../../include/AmberTopology/atom.h"
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
        double energy = cmap.return_energy_bicubic(phi, psi, resolution, topology.get_Cmap_Coefficient_Matrix_bicubic_spline(), full_grid);

        cmap.set_angles(phi, psi);
        cmap.set_energy(energy);
        cmap_energy += energy;

    }

    double LJ_energy = 0;
    std::vector<Atom>& atom_list = topology.get_atoms();
    for (size_t atomAIndex = 0; atomAIndex < topology.get_num_atoms(); ++atomAIndex)
    {
        Atom& atomA = atom_list[atomAIndex];
        for (size_t atomBIndex = atomAIndex + 1; atomBIndex < topology.get_num_atoms(); ++atomBIndex)
        {
            Atom& atomB = atom_list[atomBIndex];

            bool is14 = topology.is_14_pair(static_cast<int>(atomAIndex), static_cast<int>(atomBIndex));

            bool excluded =
                std::binary_search(atomA.excluded_atoms.begin(), atomA.excluded_atoms.end(), static_cast<int>(atomBIndex));
            // (optionally also check atomB.excluded_atoms if you don’t trust symmetry)

            if (excluded && !is14) continue;

            const double x1 = coordinates[3*atomAIndex], y1 = coordinates[3*atomAIndex + 1], z1 = coordinates[3*atomAIndex + 2];
            const double x2 = coordinates[3*atomBIndex], y2 = coordinates[3*atomBIndex + 1], z2 = coordinates[3*atomBIndex + 2];

            // double distance_AB = distance(x1, y1, z1, x2, y2, z2);
            const double dx = x1 - x2;
            const double dy = y1 - y2;
            const double dz = z1 - z2;
            double r2 = dx*dx + dy*dy + dz*dz;

            if (r2 < 1e-12) continue;

            int ti = static_cast<int>(topology.atom_list_[atomAIndex].atom_type_index) - 1;
            int tj = static_cast<int>(topology.atom_list_[atomBIndex].atom_type_index) - 1;

            std::vector<std::vector<unsigned long int>>& nbmatrix = topology.get_nb_matrix();
            unsigned long nb = nbmatrix[ti][tj];
            if (nb == 0) continue;

            double Aij = 0, Bij = 0;
            size_t p = static_cast<size_t>(nb - 1);

            if (!is14) {
                Aij = topology.get_lennard_jones_Acoefs_()[p];
                Bij = topology.get_lennard_jones_Bcoefs_()[p];
            } else {
                Aij = topology.get_lennard_jones_14_Acoefs_()[p];
                Bij = topology.get_lennard_jones_14_Bcoefs_()[p];
            }

            double energy = LennardJones::CalculateEnergy(r2, Aij, Bij);
            LJ_energy += energy;
        }
    }
    std::cout << "\n Bond:" << bond_energies << " Angle: " << angle_energy << " CosineDihedral: " << dihedral_energy <<
    " Urey-Bradley: " << urey_bradley_energy << " Impropers: " << improper_energy << " CMAP energy: " << cmap_energy << "\n"
    << "VDW: " << LJ_energy << "\n";

    calculate_forces();
    std::cout << "Calculated Forces! \n";
    std::cout << forces[0] << " " << forces[1] << " " << forces[2] << "\n";
    std::cout << forces[3] << " " << forces[4] << " " << forces[5] << "\n";
    std::cout << forces[6] << " " << forces[7] << " " << forces[8] << "\n";
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

void System::calculate_forces()
{
    clear_forces();
    calculate_forces_bonds();
    calculate_forces_UBbonds();
    calculate_forces_angles();
    calculate_forces_cosinedihedrals();
    calculate_forces_harmonicImpropers();
    calculate_forces_cmap();
    calculate_forces_LJ();
}

void System::calculate_forces_bonds() {
    const std::vector<double>& coordinates = topology.get_coordinates();

    for (auto& bond: topology.get_harmonic_bonds())
    {
        int atomAIndex = bond.get_atomA_index();
        int atomBIndex = bond.get_atomB_index();
        const double x1 = coordinates[3*atomAIndex], y1 = coordinates[3*atomAIndex + 1], z1 = coordinates[3*atomAIndex + 2];
        const double x2 = coordinates[3*atomBIndex], y2 = coordinates[3*atomBIndex + 1], z2 = coordinates[3*atomBIndex + 2];

        const double dx1 = x1 - x2;
        const double dy1 = y1 - y2;
        const double dz1 = z1 - z2;

        const double r = std::sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);

        const double k = bond.get_Bond_force_constant();
        const double r0 = bond.get_Bond_equil_length();

        if (r < 1e-12) continue;

        const double famag = -2 * k * (r - r0);
        const double fax = famag*dx1/r;
        const double fay = famag*dy1/r;
        const double faz = famag*dz1/r;

        forces[3*atomAIndex] += fax; forces[3*atomAIndex+1] += fay; forces[3*atomAIndex+2] += faz;
        forces[3*atomBIndex] -= fax; forces[3*atomBIndex+1] -= fay; forces[3*atomBIndex+2] -= faz;
    }
}

void System::calculate_forces_UBbonds() {
    // clear_forces(); // Reset gradients to 0
    const std::vector<double>& coordinates = topology.get_coordinates();

    for (auto& bond: topology.get_harmonic_UBs())
    {
        int atomAIndex = bond.get_atomA_index();
        int atomBIndex = bond.get_atomB_index();
        const double x1 = coordinates[3*atomAIndex], y1 = coordinates[3*atomAIndex + 1], z1 = coordinates[3*atomAIndex + 2];
        const double x2 = coordinates[3*atomBIndex], y2 = coordinates[3*atomBIndex + 1], z2 = coordinates[3*atomBIndex + 2];

        const double dx1 = x1 - x2;
        const double dy1 = y1 - y2;
        const double dz1 = z1 - z2;

        const double r = std::sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);

        const double k = bond.get_UB_force_constant();
        const double r0 = bond.get_UB_equil_value();

        if (r < 1e-12) continue;

        const double famag = -2 * k * (r - r0);
        const double fax = famag*dx1/r;
        const double fay = famag*dy1/r;
        const double faz = famag*dz1/r;

        forces[3*atomAIndex] += fax; forces[3*atomAIndex+1] += fay; forces[3*atomAIndex+2] += faz;
        forces[3*atomBIndex] += -fax; forces[3*atomBIndex+1] += -fay; forces[3*atomBIndex+2] += -faz;
    }
}

void System::calculate_forces_angles() {
//    clear_forces(); // Reset gradients to 0
    const std::vector<double>& coordinates = topology.get_coordinates();

    for (auto& ang: topology.get_harmonic_angles())
    {
        int atomAIndex = ang.get_atomA_index();
        int atomBIndex = ang.get_atomB_index();
        int atomCIndex = ang.get_atomC_index();

        const double x1 = coordinates[3*atomAIndex], y1 = coordinates[3*atomAIndex + 1], z1 = coordinates[3*atomAIndex + 2];
        const double x2 = coordinates[3*atomBIndex], y2 = coordinates[3*atomBIndex + 1], z2 = coordinates[3*atomBIndex + 2];
        const double x3 = coordinates[3*atomCIndex], y3 = coordinates[3*atomCIndex + 1], z3 = coordinates[3*atomCIndex + 2];

        // ba
        const double ba_x = x1 - x2;
        const double ba_y = y1 - y2;
        const double ba_z = z1 - z2;

        // bc
        const double bc_x = x3 - x2;
        const double bc_y = y3 - y2;
        const double bc_z = z3 - z2;

        const double mod_ab_sq = ba_x*ba_x + ba_y*ba_y + ba_z*ba_z;
        const double mod_ab = std::sqrt(mod_ab_sq);

        if (mod_ab < 1e-12) continue;

        // pa = norm(ba x (ba x bc))
        //    = norm ((ba.bc)ba - (modba*modba)bc)

        const double dot_ba_bc = ba_x*bc_x + ba_y*bc_y + ba_z*bc_z;
        double pa_x = dot_ba_bc*ba_x - mod_ab_sq*bc_x;
        double pa_y = dot_ba_bc*ba_y - mod_ab_sq*bc_y;
        double pa_z = dot_ba_bc*ba_z - mod_ab_sq*bc_z;

        const double mod_pa = std::sqrt(pa_x*pa_x + pa_y*pa_y + pa_z*pa_z);
        if (mod_pa < 1e-12) continue;
        pa_x /= mod_pa;
        pa_y /= mod_pa;
        pa_z /= mod_pa;

        const double mod_bc_sq = bc_x*bc_x + bc_y*bc_y + bc_z*bc_z;
        const double mod_bc = std::sqrt(mod_bc_sq);

        if (mod_bc < 1e-12) continue;

        // pc = norm(cb x (ba x bc))
        //    = norm ((bc.ba)bc -(modbc*modbc)ba)

        double pc_x = dot_ba_bc*bc_x - mod_bc_sq*ba_x;
        double pc_y = dot_ba_bc*bc_y - mod_bc_sq*ba_y;
        double pc_z = dot_ba_bc*bc_z - mod_bc_sq*ba_z;

        const double mod_pc = std::sqrt(pc_x*pc_x + pc_y*pc_y + pc_z*pc_z);
        if (mod_pc < 1e-12) continue;
        pc_x /= mod_pc;
        pc_y /= mod_pc;
        pc_z /= mod_pc;


        const double k = ang.get_Angle_force_constant();
        const double theta0 = ang.get_Angle_equil_angle();
        double theta = angle(x1, y1, z1, x2, y2, z2, x3, y3, z3);

        const double fmag = -2 * k * (theta - theta0);
        const double fax = fmag/mod_ab*pa_x;
        const double fay = fmag/mod_ab*pa_y;
        const double faz = fmag/mod_ab*pa_z;

        const double fcx = fmag/mod_bc*pc_x;
        const double fcy = fmag/mod_bc*pc_y;
        const double fcz = fmag/mod_bc*pc_z;

        forces[3*atomAIndex] += fax; forces[3*atomAIndex+1] += fay; forces[3*atomAIndex+2] += faz;
        forces[3*atomCIndex] += fcx; forces[3*atomCIndex+1] += fcy; forces[3*atomCIndex+2] += fcz;
        forces[3*atomBIndex] += -fax-fcx; forces[3*atomBIndex+1] += -fay-fcy; forces[3*atomBIndex+2] += -faz-fcz;

    }
}

void System::calculate_forces_cosinedihedrals() {
    const std::vector<double>& coordinates = topology.get_coordinates();

    for (auto& dih: topology.get_cosine_dihedrals()) {
        int atomAIndex = dih.get_atomA_index();
        int atomBIndex = dih.get_atomB_index();
        int atomCIndex = dih.get_atomC_index();
        int atomDIndex = dih.get_atomD_index();

        const double x1 = coordinates[3 * atomAIndex], y1 = coordinates[3 * atomAIndex + 1], z1 = coordinates[
                3 * atomAIndex + 2];
        const double x2 = coordinates[3 * atomBIndex], y2 = coordinates[3 * atomBIndex + 1], z2 = coordinates[
                3 * atomBIndex + 2];
        const double x3 = coordinates[3 * atomCIndex], y3 = coordinates[3 * atomCIndex + 1], z3 = coordinates[
                3 * atomCIndex + 2];
        const double x4 = coordinates[3 * atomDIndex], y4 = coordinates[3 * atomDIndex + 1], z4 = coordinates[
                3 * atomDIndex + 2];

        // ba
        const double ba_x = x1 - x2;
        const double ba_y = y1 - y2;
        const double ba_z = z1 - z2;

        // cb
        const double cb_x = x2 - x3;
        const double cb_y = y2 - y3;
        const double cb_z = z2 - z3;

        // dc
        const double dc_x = x3 - x4;
        const double dc_y = y3 - y4;
        const double dc_z = z3 - z4;


        const double k = dih.get_Dihedral_force_constant();
        const double n = dih.get_Dihedral_Periodicity();
        const double delta = dih.get_Dihedral_phase();
        double phi = dihedral(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4);
        const double tau = k * n * std::sin(n * phi - delta);

        // m = ba x cb
        const double m_x = ba_y*cb_z - ba_z*cb_y;
        const double m_y = ba_z*cb_x - ba_x*cb_z;
        const double m_z = ba_x*cb_y - ba_y*cb_x;

        // n = cb x dc

        const double n_x = cb_y*dc_z - cb_z*dc_y;
        const double n_y = cb_z*dc_x - cb_x*dc_z;
        const double n_z = cb_x*dc_y - cb_y*dc_x;

        const double mod_m_sq = m_x*m_x + m_y*m_y + m_z*m_z;
        if (mod_m_sq < 1e-24) continue;

        const double mod_n_sq = n_x*n_x + n_y*n_y + n_z*n_z;
        if (mod_n_sq < 1e-24) continue;

        const double mod_bc_sq = cb_x*cb_x + cb_y*cb_y + cb_z*cb_z;
        if (mod_bc_sq < 1e-24) continue;

        const double mod_bc = std::sqrt(mod_bc_sq);

        const double dot_ba_cb = ba_x*cb_x + ba_y*cb_y + ba_z*cb_z;
        const double dot_dc_cb = dc_x*cb_x + dc_y*cb_y + dc_z*cb_z;

        const double fax = tau*mod_bc/mod_m_sq*m_x;
        const double fay = tau*mod_bc/mod_m_sq*m_y;
        const double faz = tau*mod_bc/mod_m_sq*m_z;

        const double fdx = -tau*mod_bc/mod_n_sq*n_x;
        const double fdy = -tau*mod_bc/mod_n_sq*n_y;
        const double fdz = -tau*mod_bc/mod_n_sq*n_z;

        const double fbx = (dot_ba_cb/mod_bc_sq - 1)*fax - dot_dc_cb/mod_bc_sq*fdx;
        const double fby = (dot_ba_cb/mod_bc_sq - 1)*fay - dot_dc_cb/mod_bc_sq*fdy;
        const double fbz = (dot_ba_cb/mod_bc_sq - 1)*faz - dot_dc_cb/mod_bc_sq*fdz;

        const double fcx = (dot_dc_cb/mod_bc_sq - 1)*fdx - dot_ba_cb/mod_bc_sq*fax;
        const double fcy = (dot_dc_cb/mod_bc_sq - 1)*fdy - dot_ba_cb/mod_bc_sq*fay;
        const double fcz = (dot_dc_cb/mod_bc_sq - 1)*fdz - dot_ba_cb/mod_bc_sq*faz;

        forces[3*atomAIndex] += fax; forces[3*atomAIndex+1] += fay; forces[3*atomAIndex+2] += faz;
        forces[3*atomBIndex] += fbx; forces[3*atomBIndex+1] += fby; forces[3*atomBIndex+2] += fbz;
        forces[3*atomCIndex] += fcx; forces[3*atomCIndex+1] += fcy; forces[3*atomCIndex+2] += fcz;
        forces[3*atomDIndex] += fdx; forces[3*atomDIndex+1] += fdy; forces[3*atomDIndex+2] += fdz;
    }
}

void System::calculate_forces_harmonicImpropers() {
    const std::vector<double>& coordinates = topology.get_coordinates();

    for (auto& dih: topology.get_harmonic_impropers()) {
        int atomAIndex = dih.get_atomA_index();
        int atomBIndex = dih.get_atomB_index();
        int atomCIndex = dih.get_atomC_index();
        int atomDIndex = dih.get_atomD_index();

        const double x1 = coordinates[3 * atomAIndex], y1 = coordinates[3 * atomAIndex + 1], z1 = coordinates[
                3 * atomAIndex + 2];
        const double x2 = coordinates[3 * atomBIndex], y2 = coordinates[3 * atomBIndex + 1], z2 = coordinates[
                3 * atomBIndex + 2];
        const double x3 = coordinates[3 * atomCIndex], y3 = coordinates[3 * atomCIndex + 1], z3 = coordinates[
                3 * atomCIndex + 2];
        const double x4 = coordinates[3 * atomDIndex], y4 = coordinates[3 * atomDIndex + 1], z4 = coordinates[
                3 * atomDIndex + 2];

        // ba
        const double ba_x = x1 - x2;
        const double ba_y = y1 - y2;
        const double ba_z = z1 - z2;

        // cb
        const double cb_x = x2 - x3;
        const double cb_y = y2 - y3;
        const double cb_z = z2 - z3;

        // dc
        const double dc_x = x3 - x4;
        const double dc_y = y3 - y4;
        const double dc_z = z3 - z4;


        const double k = dih.get_IMP_force_constant();
        const double psi0 = dih.get_IMP_phase_value();
        double psi = dihedral(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4);

        double delta = psi - psi0;
        while (delta <= -M_PI) delta += 2.0 * M_PI;
        while (delta > M_PI)   delta -= 2.0 * M_PI;
        const double tau = -2.0 * k * delta;

        // m = ba x cb
        const double m_x = ba_y*cb_z - ba_z*cb_y;
        const double m_y = ba_z*cb_x - ba_x*cb_z;
        const double m_z = ba_x*cb_y - ba_y*cb_x;

        // n = cb x dc

        const double n_x = cb_y*dc_z - cb_z*dc_y;
        const double n_y = cb_z*dc_x - cb_x*dc_z;
        const double n_z = cb_x*dc_y - cb_y*dc_x;

        const double mod_m_sq = m_x*m_x + m_y*m_y + m_z*m_z;
        if (mod_m_sq < 1e-24) continue;

        const double mod_n_sq = n_x*n_x + n_y*n_y + n_z*n_z;
        if (mod_n_sq < 1e-24) continue;

        const double mod_bc_sq = cb_x*cb_x + cb_y*cb_y + cb_z*cb_z;
        if (mod_bc_sq < 1e-24) continue;

        const double mod_bc = std::sqrt(mod_bc_sq);

        const double dot_ba_cb = ba_x*cb_x + ba_y*cb_y + ba_z*cb_z;
        const double dot_dc_cb = dc_x*cb_x + dc_y*cb_y + dc_z*cb_z;

        const double fax = tau*mod_bc/mod_m_sq*m_x;
        const double fay = tau*mod_bc/mod_m_sq*m_y;
        const double faz = tau*mod_bc/mod_m_sq*m_z;

        const double fdx = -tau*mod_bc/mod_n_sq*n_x;
        const double fdy = -tau*mod_bc/mod_n_sq*n_y;
        const double fdz = -tau*mod_bc/mod_n_sq*n_z;

        const double fbx = (dot_ba_cb/mod_bc_sq - 1)*fax - dot_dc_cb/mod_bc_sq*fdx;
        const double fby = (dot_ba_cb/mod_bc_sq - 1)*fay - dot_dc_cb/mod_bc_sq*fdy;
        const double fbz = (dot_ba_cb/mod_bc_sq - 1)*faz - dot_dc_cb/mod_bc_sq*fdz;

        const double fcx = (dot_dc_cb/mod_bc_sq - 1)*fdx - dot_ba_cb/mod_bc_sq*fax;
        const double fcy = (dot_dc_cb/mod_bc_sq - 1)*fdy - dot_ba_cb/mod_bc_sq*fay;
        const double fcz = (dot_dc_cb/mod_bc_sq - 1)*fdz - dot_ba_cb/mod_bc_sq*faz;

        forces[3*atomAIndex] += fax; forces[3*atomAIndex+1] += fay; forces[3*atomAIndex+2] += faz;
        forces[3*atomBIndex] += fbx; forces[3*atomBIndex+1] += fby; forces[3*atomBIndex+2] += fbz;
        forces[3*atomCIndex] += fcx; forces[3*atomCIndex+1] += fcy; forces[3*atomCIndex+2] += fcz;
        forces[3*atomDIndex] += fdx; forces[3*atomDIndex+1] += fdy; forces[3*atomDIndex+2] += fdz;
    }
}

void System::calculate_forces_cmap()
{
    const std::vector<double>& coordinates = topology.get_coordinates();
    for (auto& cmap : topology.get_cmaps()) {

        int atomAIndex = cmap.get_atomA_index();
        int atomBIndex = cmap.get_atomB_index();
        int atomCIndex = cmap.get_atomC_index();
        int atomDIndex = cmap.get_atomD_index();
        int atomEIndex = cmap.get_atomE_index();

        const double x1 = coordinates[3 * atomAIndex], y1 = coordinates[3 * atomAIndex + 1], z1 = coordinates[
                3 * atomAIndex + 2];
        const double x2 = coordinates[3 * atomBIndex], y2 = coordinates[3 * atomBIndex + 1], z2 = coordinates[
                3 * atomBIndex + 2];
        const double x3 = coordinates[3 * atomCIndex], y3 = coordinates[3 * atomCIndex + 1], z3 = coordinates[
                3 * atomCIndex + 2];
        const double x4 = coordinates[3 * atomDIndex], y4 = coordinates[3 * atomDIndex + 1], z4 = coordinates[
                3 * atomDIndex + 2];
        const double x5 = coordinates[3 * atomEIndex], y5 = coordinates[3 * atomEIndex + 1], z5 = coordinates[
                3 * atomEIndex + 2];


        double phi = dihedral(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4);
        double psi = dihedral(x2, y2, z2, x3, y3, z3, x4, y4, z4, x5, y5, z5);

        // Calculate slopes (gradients)
        int set_index = cmap.get_parameter_set();
        const std::vector<int>& resolutions = topology.get_charmm_cmap_resolutions();
        int resolution = resolutions[set_index - 1];
        const std::vector<double>& full_grid = topology.get_cmap_grid_data(set_index);

        std::pair<double, double> grads = cmap.return_gradient_bicubic(phi, psi, resolution, topology.get_Cmap_Coefficient_Matrix_bicubic_spline(), full_grid);
        double dEdPhi = grads.first;
        double dEdPsi = grads.second;



        const double ba_x = x1 - x2;
        const double ba_y = y1 - y2;
        const double ba_z = z1 - z2;

        // cb
        const double cb_x = x2 - x3;
        const double cb_y = y2 - y3;
        const double cb_z = z2 - z3;

        // dc
        const double dc_x = x3 - x4;
        const double dc_y = y3 - y4;
        const double dc_z = z3 - z4;

        //ed
        const double ed_x = x4 - x5;
        const double ed_y = y4 - y5;
        const double ed_z = z4 - z5;


        // m = ba x cb
        const double m_x = ba_y*cb_z - ba_z*cb_y;
        const double m_y = ba_z*cb_x - ba_x*cb_z;
        const double m_z = ba_x*cb_y - ba_y*cb_x;

        // n = cb x dc

        const double n_x = cb_y*dc_z - cb_z*dc_y;
        const double n_y = cb_z*dc_x - cb_x*dc_z;
        const double n_z = cb_x*dc_y - cb_y*dc_x;

        // o = dc x ed
        const double o_x = dc_y*ed_z - dc_z*ed_y;
        const double o_y = dc_z*ed_x - dc_x*ed_z;
        const double o_z = dc_x*ed_y - dc_y*ed_x;


        const double mod_m_sq = m_x*m_x + m_y*m_y + m_z*m_z;
        if (mod_m_sq < 1e-24) continue;

        const double mod_n_sq = n_x*n_x + n_y*n_y + n_z*n_z;
        if (mod_n_sq < 1e-24) continue;

        const double mod_o_sq = o_x*o_x + o_y*o_y + o_z*o_z;
        if (mod_o_sq < 1e-24) continue;

        const double mod_bc_sq = cb_x*cb_x + cb_y*cb_y + cb_z*cb_z;
        if (mod_bc_sq < 1e-24) continue;

        const double mod_cd_sq = dc_x*dc_x + dc_y*dc_y + dc_z*dc_z;
        if (mod_cd_sq < 1e-24) continue;

        const double mod_bc = std::sqrt(mod_bc_sq);
        const double mod_cd = std::sqrt(mod_cd_sq);

        const double dot_ba_cb = ba_x*cb_x + ba_y*cb_y + ba_z*cb_z;
        const double dot_dc_cb = dc_x*cb_x + dc_y*cb_y + dc_z*cb_z;

        const double dot_cb_dc = cb_x*dc_x + cb_y*dc_y + cb_z*dc_z;
        const double dot_ed_dc = ed_x*dc_x + ed_y*dc_y + ed_z*dc_z;




        const double fax_phi = -dEdPhi*mod_bc/mod_m_sq*m_x;
        const double fay_phi = -dEdPhi*mod_bc/mod_m_sq*m_y;
        const double faz_phi = -dEdPhi*mod_bc/mod_m_sq*m_z;

        const double fdx_phi = dEdPhi*mod_bc/mod_n_sq*n_x;
        const double fdy_phi = dEdPhi*mod_bc/mod_n_sq*n_y;
        const double fdz_phi = dEdPhi*mod_bc/mod_n_sq*n_z;

        const double fbx_phi = (dot_ba_cb/mod_bc_sq - 1)*fax_phi - dot_dc_cb/mod_bc_sq*fdx_phi;
        const double fby_phi = (dot_ba_cb/mod_bc_sq - 1)*fay_phi - dot_dc_cb/mod_bc_sq*fdy_phi;
        const double fbz_phi = (dot_ba_cb/mod_bc_sq - 1)*faz_phi - dot_dc_cb/mod_bc_sq*fdz_phi;

        const double fcx_phi = (dot_dc_cb/mod_bc_sq - 1)*fdx_phi - dot_ba_cb/mod_bc_sq*fax_phi;
        const double fcy_phi = (dot_dc_cb/mod_bc_sq - 1)*fdy_phi - dot_ba_cb/mod_bc_sq*fay_phi;
        const double fcz_phi = (dot_dc_cb/mod_bc_sq - 1)*fdz_phi - dot_ba_cb/mod_bc_sq*faz_phi;


        forces[3*atomAIndex] += fax_phi; forces[3*atomAIndex+1] += fay_phi; forces[3*atomAIndex+2] += faz_phi;
        forces[3*atomBIndex] += fbx_phi; forces[3*atomBIndex+1] += fby_phi; forces[3*atomBIndex+2] += fbz_phi;
        forces[3*atomCIndex] += fcx_phi; forces[3*atomCIndex+1] += fcy_phi; forces[3*atomCIndex+2] += fcz_phi;
        forces[3*atomDIndex] += fdx_phi; forces[3*atomDIndex+1] += fdy_phi; forces[3*atomDIndex+2] += fdz_phi;

        const double fbx_psi = -dEdPsi*mod_cd/mod_n_sq*n_x;
        const double fby_psi = -dEdPsi*mod_cd/mod_n_sq*n_y;
        const double fbz_psi = -dEdPsi*mod_cd/mod_n_sq*n_z;

        const double fex_psi = dEdPsi*mod_cd/mod_o_sq*o_x;
        const double fey_psi = dEdPsi*mod_cd/mod_o_sq*o_y;
        const double fez_psi = dEdPsi*mod_cd/mod_o_sq*o_z;

        const double fcx_psi = (dot_cb_dc/mod_cd_sq - 1)*fbx_psi - dot_ed_dc/mod_cd_sq*fex_psi;
        const double fcy_psi = (dot_cb_dc/mod_cd_sq - 1)*fby_psi - dot_ed_dc/mod_cd_sq*fey_psi;
        const double fcz_psi = (dot_cb_dc/mod_cd_sq - 1)*fbz_psi - dot_ed_dc/mod_cd_sq*fez_psi;

        const double fdx_psi = (dot_ed_dc/mod_cd_sq - 1)*fex_psi - dot_cb_dc/mod_cd_sq*fbx_psi;
        const double fdy_psi = (dot_ed_dc/mod_cd_sq - 1)*fey_psi - dot_cb_dc/mod_cd_sq*fby_psi;
        const double fdz_psi = (dot_ed_dc/mod_cd_sq - 1)*fez_psi - dot_cb_dc/mod_cd_sq*fbz_psi;

        forces[3*atomBIndex] += fbx_psi; forces[3*atomBIndex+1] += fby_psi; forces[3*atomBIndex+2] += fbz_psi;
        forces[3*atomCIndex] += fcx_psi; forces[3*atomCIndex+1] += fcy_psi; forces[3*atomCIndex+2] += fcz_psi;
        forces[3*atomDIndex] += fdx_psi; forces[3*atomDIndex+1] += fdy_psi; forces[3*atomDIndex+2] += fdz_psi;
        forces[3*atomEIndex] += fex_psi; forces[3*atomEIndex+1] += fey_psi; forces[3*atomEIndex+2] += fez_psi;
    }
}

void System::calculate_forces_LJ()
{
    const std::vector<double>& coordinates = topology.get_coordinates();
    std::vector<Atom>& atom_list = topology.get_atoms();
    for (size_t atomAIndex = 0; atomAIndex < topology.get_num_atoms(); ++atomAIndex)
    {
        Atom& atomA = atom_list[atomAIndex];
        for (size_t atomBIndex = atomAIndex + 1; atomBIndex < topology.get_num_atoms(); ++atomBIndex)
        {
            Atom& atomB = atom_list[atomBIndex];

            bool is14 = topology.is_14_pair(static_cast<int>(atomAIndex), static_cast<int>(atomBIndex));

            bool excluded =
                std::binary_search(atomA.excluded_atoms.begin(), atomA.excluded_atoms.end(), static_cast<int>(atomBIndex));
            // (optionally also check atomB.excluded_atoms if you don’t trust symmetry)

            if (excluded && !is14) continue;

            const double x1 = coordinates[3*atomAIndex], y1 = coordinates[3*atomAIndex + 1], z1 = coordinates[3*atomAIndex + 2];
            const double x2 = coordinates[3*atomBIndex], y2 = coordinates[3*atomBIndex + 1], z2 = coordinates[3*atomBIndex + 2];

            // double distance_AB = distance(x1, y1, z1, x2, y2, z2);
            const double dx = x1 - x2;
            const double dy = y1 - y2;
            const double dz = z1 - z2;
            double r2 = dx*dx + dy*dy + dz*dz;

            if (r2 < 1e-12) continue;

            int ti = static_cast<int>(topology.atom_list_[atomAIndex].atom_type_index) - 1;
            int tj = static_cast<int>(topology.atom_list_[atomBIndex].atom_type_index) - 1;

            const std::vector<std::vector<unsigned long int>>& nbmatrix = topology.get_nb_matrix();
            unsigned long nb = nbmatrix[ti][tj];
            if (nb == 0) continue;

            double Aij = 0, Bij = 0;
            auto p = static_cast<size_t>(nb - 1);

            if (!is14) {
                Aij = topology.get_lennard_jones_Acoefs_()[p];
                Bij = topology.get_lennard_jones_Bcoefs_()[p];
            } else {
                Aij = topology.get_lennard_jones_14_Acoefs_()[p];
                Bij = topology.get_lennard_jones_14_Bcoefs_()[p];
            }

            double gradient = LennardJones::CalculateGradient(r2, Aij, Bij);
            const double fax = gradient * dx;
            const double fay = gradient * dy;
            const double faz = gradient * dz;

            const double fbx = -fax;
            const double fby = -fay;
            const double fbz = -faz;

            forces[3*atomAIndex] += fax; forces[3*atomAIndex+1] += fay; forces[3*atomAIndex+2] += faz;
            forces[3*atomBIndex] += fbx; forces[3*atomBIndex+1] += fby; forces[3*atomBIndex+2] += fbz;

        }
    }
}

