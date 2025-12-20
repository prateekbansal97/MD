//
// Created by Prateek Bansal on 12/20/25.
//

#include "System/System.h"
#include "System/Metrics.h"

namespace md
{
    void System::calculate_bond_energy()
    {
        const std::vector<double>& coordinates = topology_.get_coordinates();

        bond_energies_ = 0;
        for (auto& bond: topology_.get_harmonic_bonds())
        {
            const int atomAIndex = bond.get_atomA_index();
            const int atomBIndex = bond.get_atomB_index();
            const double x1 = coordinates[3*atomAIndex], y1 = coordinates[3*atomAIndex + 1], z1 = coordinates[3*atomAIndex + 2];
            const double x2 = coordinates[3*atomBIndex], y2 = coordinates[3*atomBIndex + 1], z2 = coordinates[3*atomBIndex + 2];
            const double dist = Metrics::distance(x1, y1, z1, x2, y2, z2);
            bond.set_distance(dist);

            const double energy = bond.calculate_energy(dist);
            bond.set_energy(energy);
            bond_energies_ += energy;
        }
    }

    void System::calculate_angle_energy()
    {
        const std::vector<double>& coordinates = topology_.get_coordinates();
        angle_energy_ = 0;
        for (auto& ang: topology_.get_harmonic_angles())
        {
            const int atomAIndex = ang.get_atomA_index();
            const int atomBIndex = ang.get_atomB_index();
            const int atomCIndex = ang.get_atomC_index();
            const double x1 = coordinates[3*atomAIndex], y1 = coordinates[3*atomAIndex + 1], z1 = coordinates[3*atomAIndex + 2];
            const double x2 = coordinates[3*atomBIndex], y2 = coordinates[3*atomBIndex + 1], z2 = coordinates[3*atomBIndex + 2];
            const double x3 = coordinates[3*atomCIndex], y3 = coordinates[3*atomCIndex + 1], z3 = coordinates[3*atomCIndex + 2];
            const double angle_m = Metrics::angle(x1, y1, z1, x2, y2, z2, x3, y3, z3);
            ang.set_angle(angle_m);

            const double energy = ang.calculate_energy(angle_m);
            ang.set_energy(energy);
            angle_energy_ += energy;
        }
    }

    void System::calculate_dihedral_energy()
    {
        const std::vector<double>& coordinates = topology_.get_coordinates();

        dihedral_energy_ = 0;
        for (auto& dih: topology_.get_cosine_dihedrals())
        {
            const int atomAIndex = dih.get_atomA_index();
            const int atomBIndex = dih.get_atomB_index();
            const int atomCIndex = dih.get_atomC_index();
            const int atomDIndex = dih.get_atomD_index();
            const double x1 = coordinates[3*atomAIndex], y1 = coordinates[3*atomAIndex + 1], z1 = coordinates[3*atomAIndex + 2];
            const double x2 = coordinates[3*atomBIndex], y2 = coordinates[3*atomBIndex + 1], z2 = coordinates[3*atomBIndex + 2];
            const double x3 = coordinates[3*atomCIndex], y3 = coordinates[3*atomCIndex + 1], z3 = coordinates[3*atomCIndex + 2];
            const double x4 = coordinates[3*atomDIndex], y4 = coordinates[3*atomDIndex + 1], z4 = coordinates[3*atomDIndex + 2];
            const double angle_m = Metrics::dihedral(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4);
            dih.set_cosine_dihedral(angle_m);
            const double energy = dih.calculate_energy(angle_m);
            dih.set_energy(energy);
            dihedral_energy_ += energy;
        }
    }

    void System::calculate_UB_energy()
    {
        const std::vector<double>& coordinates = topology_.get_coordinates();

        urey_bradley_energy_ = 0;
        for (auto& ub : topology_.get_harmonic_UBs())
        {
            const int a1 = ub.get_atomA_index();
            const int a2 = ub.get_atomB_index();


            const double dist = Metrics::distance(
                    coordinates[3*a1], coordinates[3*a1+1], coordinates[3*a1+2],
                    coordinates[3*a2], coordinates[3*a2+1], coordinates[3*a2+2]
            );

            ub.set_distance_value(dist);
            const double energy = ub.calculate_energy(dist);
            ub.set_energy(energy);
            urey_bradley_energy_ += energy;
        }
    }

    void System::calculate_improper_energy()
    {
        const std::vector<double>& coordinates = topology_.get_coordinates();

        improper_energy_ = 0;
        for (auto& imp : topology_.get_harmonic_impropers())
        {
            const int a1 = imp.get_atomA_index();
            const int a2 = imp.get_atomB_index();
            const int a3 = imp.get_atomC_index();
            const int a4 = imp.get_atomD_index();

            const double angle_m = Metrics::dihedral(
                    coordinates[3*a1], coordinates[3*a1+1], coordinates[3*a1+2],
                    coordinates[3*a2], coordinates[3*a2+1], coordinates[3*a2+2],
                    coordinates[3*a3], coordinates[3*a3+1], coordinates[3*a3+2],
                    coordinates[3*a4], coordinates[3*a4+1], coordinates[3*a4+2]
            );

            // Assuming you added a 'current_angle' member or setter to HarmonicImproper
            imp.set_imp_dihedral(angle_m);
            const double energy = imp.calculate_energy(angle_m);
            imp.set_energy(energy);
            improper_energy_ += energy;
        }
    }

    void System::calculate_CMAP_energy()
    {
        const std::vector<double>& coordinates = topology_.get_coordinates();

        CMAP_energy_ = 0.0;
        for (auto& cmap : topology_.get_cmaps())
        {
            const int a1 = cmap.get_atomA_index();
            const int a2 = cmap.get_atomB_index();
            const int a3 = cmap.get_atomC_index();
            const int a4 = cmap.get_atomD_index();
            const int a5 = cmap.get_atomE_index();

            // Calculate First Torsion (A-B-C-D)
            const double phi = Metrics::dihedral(
                    coordinates[3*a1], coordinates[3*a1+1], coordinates[3*a1+2],
                    coordinates[3*a2], coordinates[3*a2+1], coordinates[3*a2+2],
                    coordinates[3*a3], coordinates[3*a3+1], coordinates[3*a3+2],
                    coordinates[3*a4], coordinates[3*a4+1], coordinates[3*a4+2]
            );

            // Calculate Second Torsion (B-C-D-E)
            const double psi = Metrics::dihedral(
                    coordinates[3*a2], coordinates[3*a2+1], coordinates[3*a2+2],
                    coordinates[3*a3], coordinates[3*a3+1], coordinates[3*a3+2],
                    coordinates[3*a4], coordinates[3*a4+1], coordinates[3*a4+2],
                    coordinates[3*a5], coordinates[3*a5+1], coordinates[3*a5+2]
            );

            const int set_index = cmap.get_parameter_set();
            const std::vector<int>& resolutions = topology_.get_charmm_cmap_resolutions();
            const int resolution = resolutions[set_index - 1];
            const std::vector<double>& full_grid = topology_.get_cmap_grid_data(set_index);


            // Pass the OFFSET and the FULL GRID to your return_energy function
            const double energy = cmap.return_energy_bicubic(phi, psi, resolution, topology_.get_Cmap_Coefficient_Matrix_bicubic_spline(), full_grid);

            cmap.set_angles(phi, psi);
            cmap.set_energy(energy);
            CMAP_energy_ += energy;
        }
    }
}