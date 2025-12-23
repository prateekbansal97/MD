//
// Created by Prateek Bansal on 12/20/25.
//

#include "System/System.h"
#include "NonBonded/LennardJones.h"
#include "NonBonded/CoulombicEE.h"


namespace md
{
    void System::calculate_forces_LJ()
    {
        const std::vector<double>& coordinates = topology_.get_coordinates();
        std::vector<Atom>& atom_list = topology_.get_atom_list();

        const auto& A    = topology_.get_lennard_jones_Acoefs_();
        const auto& B    = topology_.get_lennard_jones_Bcoefs_();
        const auto& A14  = topology_.get_lennard_jones_14_Acoefs_();
        const auto& B14  = topology_.get_lennard_jones_14_Bcoefs_();


        for (size_t atomAIndex = 0; atomAIndex < natoms_; ++atomAIndex)
        {
            Atom& atomA = atom_list[atomAIndex];
            const double x1 = coordinates[3*atomAIndex], y1 = coordinates[3*atomAIndex + 1], z1 = coordinates[3*atomAIndex + 2];
            const std::vector<int>& excl = atomA.get_excluded_atoms();

            for (size_t atomBIndex = atomAIndex + 1; atomBIndex < natoms_; ++atomBIndex)
            {

                bool is14 = false;

                if (std::ranges::binary_search(excl, static_cast<int>(atomBIndex))) {
                    is14 = topology_.is_14_pair(static_cast<int>(atomAIndex), static_cast<int>(atomBIndex));
                    if (!is14) continue; // excluded and not 1-4 => skip
                }


                const double x2 = coordinates[3*atomBIndex], y2 = coordinates[3*atomBIndex + 1], z2 = coordinates[3*atomBIndex + 2];

                double dx = x1 - x2;
                double dy = y1 - y2;
                double dz = z1 - z2;
                apply_min_image(dx,dy,dz);
                const double r2 = dx*dx + dy*dy + dz*dz;

                if (r2 > lj_cutoff2_) continue;
                if (r2 < 1e-12) continue;

                const double r = std::sqrt(r2);
                if (r < 1e-12) continue;


                const int ti = type_[atomAIndex];
                const int tj = type_[atomBIndex];
                const unsigned long nb = nb_flat_[static_cast<size_t>(ti) * nTypes_ + tj];

                if (nb == 0) continue;

                double Aij = 0, Bij = 0;
                const auto p = static_cast<size_t>(nb - 1);

                if (!is14) {
                    Aij = A[p];
                    Bij = B[p];
                } else {
                    Aij = A14[p];
                    Bij = B14[p];
                }


                double gradient = NonBonded::LennardJones::CalculateGradient(r2, Aij, Bij);
                const double gcut = (!is14 ? lj_Gcut_[p] : lj_Gcut14_[p]);
                gradient -= gcut * (lj_cutoff_ / r);

                const double fax = gradient * dx;
                const double fay = gradient * dy;
                const double faz = gradient * dz;

                const double fbx = -fax;
                const double fby = -fay;
                const double fbz = -faz;

                forces_[3*atomAIndex] += fax; forces_[3*atomAIndex+1] += fay; forces_[3*atomAIndex+2] += faz;
                forces_[3*atomBIndex] += fbx; forces_[3*atomBIndex+1] += fby; forces_[3*atomBIndex+2] += fbz;
            }
        }
    }
    

    void System::calculate_forces_LJ_pairlist()
    {
        const std::vector<double>& coordinates = topology_.get_coordinates();

        const auto& A    = topology_.get_lennard_jones_Acoefs_();
        const auto& B    = topology_.get_lennard_jones_Bcoefs_();
        const auto& A14  = topology_.get_lennard_jones_14_Acoefs_();
        const auto& B14  = topology_.get_lennard_jones_14_Bcoefs_();

        for (auto& [atomAIndex, atomBIndex, p, is14]: lj_pairs_)
        {
            const double x1 = coordinates[3*atomAIndex], y1 = coordinates[3*atomAIndex + 1], z1 = coordinates[3*atomAIndex + 2];
            const double x2 = coordinates[3*atomBIndex], y2 = coordinates[3*atomBIndex + 1], z2 = coordinates[3*atomBIndex + 2];

            double dx = x1 - x2;
            double dy = y1 - y2;
            double dz = z1 - z2;
            apply_min_image(dx,dy,dz);
            const double r2 = dx*dx + dy*dy + dz*dz;

            if (r2 > lj_cutoff2_) continue;
            if (r2 < 1e-12) continue;

            const double r = std::sqrt(r2);
            if (r < 1e-12) continue;

            double Aij = 0, Bij = 0;

            if (!is14) {
                Aij = A[p];
                Bij = B[p];
            } else {
                Aij = A14[p];
                Bij = B14[p];
            }


            double gradient = NonBonded::LennardJones::CalculateGradient(r2, Aij, Bij);


            if (lj_use_switch_) {
                double S, dSdr;
                lj_switch_factors(r, S, dSdr);
                const double U = NonBonded::LennardJones::CalculateEnergy(r2, Aij, Bij);
                gradient = S*gradient - (U / r) * dSdr;
            } else {
                // Force-shift to make force go to 0 at cutoff
                const double gcut = (!is14 ? lj_Gcut_[p] : lj_Gcut14_[p]);   // g(rc) = -U'(rc)/rc
                gradient -= gcut * (lj_cutoff_ / r);                                // => (-U'(r)+U'(rc))/r
            }

            const double fax = gradient * dx;
            const double fay = gradient * dy;
            const double faz = gradient * dz;

            const double fbx = -fax;
            const double fby = -fay;
            const double fbz = -faz;

            forces_[3*atomAIndex] += fax; forces_[3*atomAIndex+1] += fay; forces_[3*atomAIndex+2] += faz;
            forces_[3*atomBIndex] += fbx; forces_[3*atomBIndex+1] += fby; forces_[3*atomBIndex+2] += fbz;
        }
    }


    void System::calculate_forces_EE()
    {
        const std::vector<double>& coordinates = topology_.get_coordinates();
        std::vector<Atom>& atom_list = topology_.get_atom_list();


        for (size_t atomAIndex = 0; atomAIndex < natoms_; ++atomAIndex)
        {
            Atom& atomA = atom_list[atomAIndex];
            const std::vector<int>& excl = atomA.get_excluded_atoms();
            const double x1 = coordinates[3*atomAIndex], y1 = coordinates[3*atomAIndex + 1], z1 = coordinates[3*atomAIndex + 2];

            for (size_t atomBIndex = atomAIndex + 1; atomBIndex < natoms_; ++atomBIndex)
            {
                bool is14 = false;

                if (std::ranges::binary_search(excl, static_cast<int>(atomBIndex))) {
                    is14 = topology_.is_14_pair(static_cast<int>(atomAIndex), static_cast<int>(atomBIndex));
                    if (!is14) continue; // excluded and not 1-4 => skip
                }

                const double x2 = coordinates[3*atomBIndex], y2 = coordinates[3*atomBIndex + 1], z2 = coordinates[3*atomBIndex + 2];

                double dx = x1 - x2;
                double dy = y1 - y2;
                double dz = z1 - z2;
                apply_min_image(dx,dy,dz);
                const double r2 = dx*dx + dy*dy + dz*dz;
                if (r2 < 1e-12) continue;
                if (r2 > lj_cutoff2_) continue;
                const double r = std::sqrt(r2);

                double scale = 1.0;
                if (is14) scale = topology_.get_scee_scale_for_pair(static_cast<int>(atomAIndex), static_cast<int>(atomBIndex));

                const double gradient = (1/scale) * NonBonded::CoulombicEE::CalculateGradient(r, q_[atomAIndex], q_[atomBIndex], 1);
                const double fax = gradient * dx;
                const double fay = gradient * dy;
                const double faz = gradient * dz;

                const double fbx = -fax;
                const double fby = -fay;
                const double fbz = -faz;

                forces_[3*atomAIndex] += fax; forces_[3*atomAIndex+1] += fay; forces_[3*atomAIndex+2] += faz;
                forces_[3*atomBIndex] += fbx; forces_[3*atomBIndex+1] += fby; forces_[3*atomBIndex+2] += fbz;
            }
        }
    }

    void System::calculate_forces_EE_pairlist()
    {
        const std::vector<double>& coordinates = topology_.get_coordinates();

        for (auto& [atomAIndex, atomBIndex, p, is14]: ee_pairs_)
        {
            const double x1 = coordinates[3*atomAIndex], y1 = coordinates[3*atomAIndex + 1], z1 = coordinates[3*atomAIndex + 2];
            const double x2 = coordinates[3*atomBIndex], y2 = coordinates[3*atomBIndex + 1], z2 = coordinates[3*atomBIndex + 2];

            double dx = x1 - x2;
            double dy = y1 - y2;
            double dz = z1 - z2;
            apply_min_image(dx,dy,dz);
            const double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 > ee_cutoff2_) continue;
            if (r2 < 1e-12) continue;

            const double r = std::sqrt(r2);
            if (r < 1e-12) continue;

            double scale = 1.0;
            if (is14) scale = topology_.get_scee_scale_for_pair(atomAIndex, atomBIndex);

            const double gradient = (1/scale) * NonBonded::CoulombicEE::CalculateGradient(r, q_[atomAIndex], q_[atomBIndex], 1);
            const double fax = gradient * dx;
            const double fay = gradient * dy;
            const double faz = gradient * dz;

            const double fbx = -fax;
            const double fby = -fay;
            const double fbz = -faz;

            forces_[3*atomAIndex] += fax; forces_[3*atomAIndex+1] += fay; forces_[3*atomAIndex+2] += faz;
            forces_[3*atomBIndex] += fbx; forces_[3*atomBIndex+1] += fby; forces_[3*atomBIndex+2] += fbz;
        }
    }
}