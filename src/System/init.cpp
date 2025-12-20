//
// Created by Prateek Bansal on 12/20/25.
//

#include "System/System.h"

namespace md
{
    void System::init()
    {

        build_nonbonded_cache();

        const auto boxdim = topology_.get_box_dimensions();
        set_box(boxdim[0], boxdim[1], boxdim[2]);

        set_lj_cutoff(10.0);
        set_lj_skin_(2.0);
        set_lj_switch_radius(9.0);

        set_ee_cutoff(10.0);
        set_ee_skin(2.0);



        setup_neighbor_rebuild_threshold_lj();
        setup_neighbor_rebuild_threshold_ee();
        init_box();

        if (const double halfmin = 0.5 * std::min({boxLx_, boxLy_, boxLz_}); pbc_enabled_ && lj_cutoff_ > halfmin) {
            std::cerr << "LJ cutoff too large for minimum image.\n";
        }

        precompute_lj_energy_force_shift();
        setup_pairlists();
        copy_ref_coords_lj();
        copy_ref_coords_ee();

        calculate_energies();
        std::cout << "\nBond:" << bond_energies_ << "\nAngle: " << angle_energy_ << "\nCosineDihedral: " << dihedral_energy_ <<
        "\nUrey-Bradley: " << urey_bradley_energy_ << "\nImpropers: " << improper_energy_ << "\nCMAP energy: " << CMAP_energy_
        << "\nVDW: " << LJ_energy_pairlist_ << "\nEE: " << EE_energy_ << "\nEE pairlist: " << EE_energy_pairlist_ << "\n";

        calculate_forces();
        std::cout << "Calculated Forces! \n"
        << forces_[0] << " " << forces_[1] << " " << forces_[2] << "\n"
        << forces_[3] << " " << forces_[4] << " " << forces_[5] << "\n"
        << forces_[6] << " " << forces_[7] << " " << forces_[8] << "\n";

    }
}
