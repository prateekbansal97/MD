//
// Created by Prateek Bansal on 12/21/25.
//

#include "AmberTopology/AmberTopology.h"

namespace md::AmberTopology
{
    void AmberTopology::process_charmm_urey_bradley_assign()
    {
        for (int i = 0; i + 3 < static_cast<int>(charmm_urey_bradley_indices_.size()); i=i+3)
        {
            const int indexA = charmm_urey_bradley_indices_[i] - 1;
            const int indexB = charmm_urey_bradley_indices_[i+1] - 1;
            const int index_CUB_type = charmm_urey_bradley_indices_[i+2] - 1;


            check_if_valid_indices("atom_list_", atom_list_, indexA, indexB);
            Bonded::HarmonicUB UB_bond = Bonded::HarmonicUB(indexA, indexB);
            UB_bond.set_type(index_CUB_type);
            HarmonicUB_list_.push_back(UB_bond);
        }
    }

    void AmberTopology::process_charmm_urey_bradley_assign_ffparams()
    {
        for (auto& UB_bond: HarmonicUB_list_)
        {
            const int type = UB_bond.get_type();
            check_if_valid_indices("charmm_urey_bradley_force_constants_", charmm_urey_bradley_force_constants_, type);
            check_if_valid_indices("charmm_urey_bradley_equil_values_", charmm_urey_bradley_equil_values_, type);
            const double force_constant = charmm_urey_bradley_force_constants_[type];
            const double equil_value = charmm_urey_bradley_equil_values_[type];
            UB_bond.set_UB_force_constant(force_constant);
            UB_bond.set_UB_equil_value(equil_value);
        }
    }

    void AmberTopology::process_charmm_improper_assign()
    {
        for (int i = 0; i + 4 < static_cast<int>(charmm_improper_indices_.size()); i=i+5)
        {
            const int indexA = charmm_improper_indices_[i] - 1;
            const int indexB = charmm_improper_indices_[i+1] - 1;
            const int indexC = charmm_improper_indices_[i+2] - 1;
            const int indexD = charmm_improper_indices_[i+3] - 1;
            const int index_CI_type = charmm_improper_indices_[i+4] - 1;

            check_if_valid_indices("atomA_index", atom_list_, indexA, indexB, indexC, indexD);
            Bonded::HarmonicImproper IMP_bond = Bonded::HarmonicImproper(indexA, indexB, indexC, indexD);
            IMP_bond.set_type(index_CI_type);
            HarmonicIMP_list_.push_back(IMP_bond);
        }
    }

    void AmberTopology::process_improper_bradley_assign_ffparams()
    {
        for (Bonded::HarmonicImproper& IMP_bond: HarmonicIMP_list_)
        {
            const int type = IMP_bond.get_type();
            check_if_valid_indices("charmm_improper_force_constants_", charmm_improper_force_constants_, type);
            check_if_valid_indices("charmm_improper_phases_", charmm_improper_phases_, type);

            const double force_constant = charmm_improper_force_constants_[type];
            const double phase_value = charmm_improper_phases_[type];

            IMP_bond.set_IMP_force_constant(force_constant);
            IMP_bond.set_IMP_phase_value(phase_value);
        }
    }

    void AmberTopology::create_bonds_including_H()
    {
        for (size_t i = 0; i + 3 < bonds_including_h_.size(); i = i + 3)
        {
            const int indexA = (bonds_including_h_[i])/3;
            const int indexB = (bonds_including_h_[i+1])/3;
            const int force_type = bonds_including_h_[i+2] - 1;

            check_if_valid_indices("atomA_index", atom_list_, indexA, indexB);

            if (force_type < 0)
            {
                std::cout << "Invalid indices";
            }

            check_if_valid_indices("bond_force_constants_", bond_force_constants_, force_type);
            check_if_valid_indices("bond_force_equils_", bond_force_equils_, force_type);

            const double force_constant = bond_force_constants_[force_type];
            const double equil_length = bond_force_equils_[force_type];

            Bonded::HarmonicBond Bond = Bonded::HarmonicBond(force_constant, equil_length, indexA, indexB, true);
            HarmonicBond_list_.push_back(Bond);
        }
    }

    void AmberTopology::create_bonds_without_H()
    {
        for (size_t i = 0; i + 3 < bonds_without_h_.size(); i = i + 3)
        {
            const int indexA = (bonds_without_h_[i])/3;
            const int indexB = (bonds_without_h_[i+1])/3;
            const int force_type = bonds_without_h_[i+2] - 1;


            check_if_valid_indices("atom_list_", atom_list_, indexA, indexB);
            if (force_type < 0)
            {
                std::cout << "Invalid indices";
            }

            check_if_valid_indices("bond_force_constants_", bond_force_constants_, force_type);
            check_if_valid_indices("bond_force_equils_", bond_force_equils_, force_type);

            const double force_constant = bond_force_constants_[force_type];
            const double equil_length = bond_force_equils_[force_type];

            Bonded::HarmonicBond Bond = Bonded::HarmonicBond(force_constant, equil_length, indexA, indexB, false);
            HarmonicBond_list_.push_back(Bond);
        }
    }

    void AmberTopology::create_angles_including_H()
    {
        for (size_t i = 0; i + 4 < angles_including_h_.size(); i = i + 4)
        {
            const int indexA = (angles_including_h_[i])/3;
            const int indexB = (angles_including_h_[i+1])/3;
            const int indexC = (angles_including_h_[i+2])/3;
            const int force_type = angles_including_h_[i+3] - 1;

            check_if_valid_indices("atom_list_", atom_list_, indexA, indexB, indexC);

            if (force_type < 0)
            {
                std::cout << "Invalid indices";
            }

            check_if_valid_indices("angle_force_constants_", angle_force_constants_, force_type);
            check_if_valid_indices("angle_force_equils_", angle_force_equils_, force_type);
            const double force_constant = angle_force_constants_[force_type];
            const double equil_angle = angle_force_equils_[force_type];

            Bonded::HarmonicAngle Angle = Bonded::HarmonicAngle(force_constant, equil_angle, indexA, indexB, indexC, true);
            HarmonicAngle_list_.push_back(Angle);
        }
    }

    void AmberTopology::create_angles_without_H()
    {
        for (size_t i = 0; i + 4 < angles_without_h_.size(); i = i + 4)
        {
            const int indexA = (angles_without_h_[i])/3;
            const int indexB = (angles_without_h_[i+1])/3;
            const int indexC = (angles_without_h_[i+2])/3;
            const int force_type = angles_without_h_[i+3] - 1;
            check_if_valid_indices("atom_list_", atom_list_, indexA, indexB, indexC);

            if (force_type < 0)
            {
                std::cout << "Invalid indices";
            }


            check_if_valid_indices("angle_force_constants_", angle_force_constants_, force_type);
            check_if_valid_indices("angle_force_equils_", angle_force_equils_, force_type);
            const double force_constant = angle_force_constants_[force_type];
            const double equil_angle = angle_force_equils_[force_type];

            Bonded::HarmonicAngle Angle = Bonded::HarmonicAngle(force_constant, equil_angle, indexA, indexB, indexC, false);
            HarmonicAngle_list_.push_back(Angle);
        }
    }

    void AmberTopology::create_dihedrals_including_H()
    {
        for (size_t i = 0; i + 5 < dihedrals_including_h_.size(); i = i + 5)
        {

            bool exclude_14 = false;
            bool improper = false;
            const int indexA = (dihedrals_including_h_[i])/3;
            const int indexB = (dihedrals_including_h_[i+1])/3;
            int indexC = (dihedrals_including_h_[i+2])/3;

            if (indexC < 0)
            {
                exclude_14 = true;
                indexC = std::abs(indexC);
            }
            int indexD = (dihedrals_including_h_[i+3])/3;
            if (indexD < 0)
            {
                improper = true;
                indexD = std::abs(indexD);
            }

            check_if_valid_indices("atom_list_", atom_list_, indexA, indexB, indexC, indexD);

            const int force_type = dihedrals_including_h_[i+4] - 1;


            check_if_valid_indices("dihedral_force_constants_", dihedral_force_constants_, force_type);
            check_if_valid_indices("dihedral_phases_", dihedral_phases_, force_type);
            check_if_valid_indices("dihedral_periodicities_", dihedral_periodicities_, force_type);

            const double force_constant = dihedral_force_constants_[force_type];
            const double phase = dihedral_phases_[force_type];
            const double periodicity = dihedral_periodicities_[force_type];

            Bonded::CosineDihedral Dihedral = Bonded::CosineDihedral(force_constant, phase, periodicity, indexA, indexB, indexC, indexD, true, improper, exclude_14);
            Dihedral.set_type(force_type);
            CosineDihedral_list_.push_back(Dihedral);
        }
    }

    void AmberTopology::create_dihedrals_without_H()
    {
        for (size_t i = 0; i + 5 < dihedrals_without_h_.size(); i = i + 5)
        {
            bool exclude_14 = false;
            bool improper = false;
            const int indexA = (dihedrals_without_h_[i])/3;
            const int indexB = (dihedrals_without_h_[i+1])/3;
            int indexC = (dihedrals_without_h_[i+2])/3;
            if (indexC < 0)
            {
                exclude_14 = true;
                indexC = std::abs(indexC);
            }
            int indexD = (dihedrals_without_h_[i+3])/3;
            if (indexD < 0)
            {
                improper = true;
                indexD = std::abs(indexD);
            }
            const int force_type = dihedrals_without_h_[i+4] - 1;

            if (force_type < 0)
            {
                std::cout << "Invalid indices";
            }

            check_if_valid_indices("atom_list_", atom_list_, indexA, indexB, indexC, indexD);
            check_if_valid_indices("dihedral_force_constants_", dihedral_force_constants_, force_type);
            check_if_valid_indices("dihedral_phases_", dihedral_phases_, force_type);
            check_if_valid_indices("dihedral_periodicities_", dihedral_periodicities_, force_type);

            const double force_constant = dihedral_force_constants_[force_type];
            const double phase = dihedral_phases_[force_type];
            const double periodicity = dihedral_periodicities_[force_type];

            Bonded::CosineDihedral Dihedral = Bonded::CosineDihedral(force_constant, phase, periodicity, indexA, indexB, indexC, indexD, false, improper, exclude_14);
            Dihedral.set_type(force_type);
            CosineDihedral_list_.push_back(Dihedral);
        }
    }

    void AmberTopology::create_Charmm_Cmap_Index_assign()
    {
        for (size_t i = 0; i + 6 < charmm_cmap_index_.size(); i = i + 6)
        {
            const int indexA = charmm_cmap_index_[i] - 1;
            const int indexB = charmm_cmap_index_[i+1] - 1;
            const int indexC = charmm_cmap_index_[i+2] - 1;
            const int indexD = charmm_cmap_index_[i+3] - 1;
            const int indexE = charmm_cmap_index_[i+4] - 1;
            check_if_valid_indices("atom_list_Charmm_Cmap_Index_assign", atom_list_, indexA, indexB, indexC, indexD, indexE);

            const int parameter_set = charmm_cmap_index_[i+5];


            Bonded::CMapGroup group = Bonded::CMapGroup(parameter_set, indexA, indexB, indexC, indexD, indexE);
            CMapGroup_list_.push_back(group);
        }
    }

}
