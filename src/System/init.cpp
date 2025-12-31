//
// Created by Prateek Bansal on 12/20/25.
//

#include "System/System.h"
// #include <fftw3.h>

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

        const double comx = calculate_center_of_mass_along_direction("x");
        const double comy = calculate_center_of_mass_along_direction("y");
        const double comz = calculate_center_of_mass_along_direction("z");
        std::cout << "Center of mass along x: " << comx << "\n" <<
                    "Center of mass along y: " << comy << "\n" <<
                    "Center of mass along z: " << comz << std::endl;

        std::abs(comx - boxLx_/2) < 1 ? std::cout << "Not moving box along x.\n" : std::cout << "Moving box along x by " << comx - boxLx_/2 << std::endl; move_box("x", comx - boxLx_/2);
        std::abs(comy - boxLy_/2) < 1 ? std::cout << "Not moving box along y.\n" : std::cout << "Moving box along y by " << comy - boxLy_/2 << std::endl; move_box("y", comy - boxLy_/2);
        std::abs(comz - boxLz_/2) < 1 ? std::cout << "Not moving box along z.\n" : std::cout << "Moving box along z by " << comz - boxLz_/2 << std::endl; move_box("z", comz - boxLz_/2);



        set_ewald_error_tolerance(1e-6);
        calculate_ewald_alpha();
        calculate_ewald_nnodes();
        calculate_pme_grid_point_positions();
        std::cout << "nX: " << newald_mesh_x << " nY: " << newald_mesh_y << " nZ: " << newald_mesh_z << std::endl;

        spread_charge_pme();
        init_pme_resources();
        perform_forward_fft();
        solve_poisson_kspace();
        perform_backward_fft();
        gather_forces_pme();

        calculate_energies();
        std::cout << "\nBond:" << bond_energies_ << "\nAngle: " << angle_energy_ << "\nCosineDihedral: " << dihedral_energy_ <<
        "\nUrey-Bradley: " << urey_bradley_energy_ << "\nImpropers: " << improper_energy_ << "\nCMAP energy: " << CMAP_energy_
        << "\nVDW: " << LJ_energy_pairlist_ << "\nEE: " << EE_energy_ << "\nEE pairlist: " << EE_energy_pairlist_ << "\n" << " EE_energy_pairlist_cutoff_"
        << EE_energy_pairlist_cutoff_ << "\n" <<
            "EE_energy_pairlist_ewald_direct_ " << EE_energy_pairlist_ewald_direct_  << "\n" <<
            "EE_energy_pairlist_ewald_self_" << EE_energy_pairlist_ewald_self_ << "\n";


        calculate_forces();
        std::cout << "Calculated Forces! \n"
        << forces_[0] << " " << forces_[1] << " " << forces_[2] << "\n"
        << forces_[3] << " " << forces_[4] << " " << forces_[5] << "\n"
        << forces_[6] << " " << forces_[7] << " " << forces_[8] << "\n";

    }
}
