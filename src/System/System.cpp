//
// Created by Prateek Bansal on 12/9/25.
//

#include "System/System.h"

namespace md
{
    void System::calculate_energies()
    {
        calculate_bond_energy();
        calculate_angle_energy();
        calculate_dihedral_energy();
        calculate_UB_energy();
        calculate_improper_energy();
        calculate_CMAP_energy();
        calculate_LJ_energy_pairlist();
        calculate_EE_energy();
        calculate_EE_energy_pairlist();
        calculate_EE_energy_pairlist_with_cutoff();
        calculate_EE_ewald_direct_term();
        calculate_EE_ewald_self_term();
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
        // calculate_forces_LJ();
        calculate_forces_LJ_pairlist();
        // calculate_forces_EE();
        calculate_forces_EE_pairlist();
    }
}// namespace md




