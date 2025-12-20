//
// Created by Prateek Bansal on 12/20/25.
//

#include "System/System.h"
#include "NonBonded/LennardJones.h"

namespace md
{

    void System::set_ee_cutoff(const double cutoff)
    {
        ee_cutoff_ = cutoff;
        ee_cutoff2_ = cutoff*cutoff;
    }

    void System::set_lj_list_cutoff(const double cutoff)
    {
        lj_list_cutoff_ = cutoff;
        lj_list_cutoff2_ = cutoff*cutoff;
    }

    void System::set_lj_skin_(const double skin)
    {
        lj_skin_ = skin;
        lj_list_cutoff_ = lj_skin_ + lj_cutoff_;
        lj_list_cutoff2_ = lj_list_cutoff_*lj_list_cutoff_;
    }

    void System::set_lj_switch_radius(const double rswitch) {
        lj_switch_ = rswitch;
        lj_use_switch_ = (rswitch > 0.0);
        // optional safety:
        if (lj_use_switch_ && !(lj_switch_ < lj_cutoff_)) {
            throw std::runtime_error("lj_switch must be < lj_cutoff");
        }
    }

    void System::set_lj_switch(const double rswitch) {
        lj_switch_ = rswitch;
        lj_use_switch_ = (rswitch > 0.0);
        // optional safety:
        if (lj_use_switch_ && !(lj_switch_ < lj_cutoff_)) {
            throw std::runtime_error("lj_switch must be < lj_cutoff");
        }
    }


    void System::setup_pairlists()
    {
        build_lj_pairlist(); // fills lj_pairs_

        if ((lj_skin_ == ee_skin_) && (lj_cutoff_ == ee_cutoff_)) {
            std::cout << "Using the same lists for LJ and EE...\n";
            ee_pairs_ = std::span<const NBPair>(lj_pairs_);
            ee_pairs_buf_.clear();              // no duplicate allocation in use
            ee_pairs_buf_.shrink_to_fit();    // optional if you want to release memory

        } else {
            build_ee_pairlist();                // fills ee_pairs_buf_
            ee_pairs_ = std::span<const NBPair>(ee_pairs_buf_);
        }
    }

    void System::setup_neighbor_rebuild_threshold_lj() {
        const double half_skin = 0.5 * lj_skin_;
        rebuild_disp2_lj_ = half_skin * half_skin;
        coordinates_ref_wrapped_lj_.resize(3 * natoms_);
    }

    void System::setup_neighbor_rebuild_threshold_ee() {
        const double half_skin = 0.5 * ee_skin_;
        rebuild_disp2_ee_ = half_skin * half_skin;
        coordinates_ref_wrapped_ee_.resize(3 * natoms_);
    }

    void System::copy_ref_coords_lj() {
        const auto& c = topology_.get_coordinates();
        for (size_t i = 0; i < 3 * natoms_; ++i) coordinates_ref_wrapped_lj_[i] = c[i];
    }

    void System::copy_ref_coords_ee() {
        const auto& c = topology_.get_coordinates();
        for (size_t i = 0; i < 3 * natoms_; ++i) coordinates_ref_wrapped_ee_[i] = c[i];
    }

    bool System::ljpairs_needs_rebuild() {
        const auto& c = topology_.get_coordinates();
        double max_dr2 = 0.0;

        for (int i = 0; i < natoms_; ++i) {
            double dx = c[3*i]   - coordinates_ref_wrapped_lj_[3*i];
            double dy = c[3*i+1] - coordinates_ref_wrapped_lj_[3*i+1];
            double dz = c[3*i+2] - coordinates_ref_wrapped_lj_[3*i+2];
            apply_min_image(dx, dy, dz);
            if (const double dr2 = dx*dx + dy*dy + dz*dz; dr2 > max_dr2) max_dr2 = dr2;
            if (max_dr2 > rebuild_disp2_lj_) return true; // early exit
        }
        return false;
    }

    bool System::eepairs_needs_rebuild() {
        const auto& c = topology_.get_coordinates();
        double max_dr2 = 0.0;

        for (int i = 0; i < natoms_; ++i) {
            double dx = c[3*i]   - coordinates_ref_wrapped_ee_[3*i];
            double dy = c[3*i+1] - coordinates_ref_wrapped_ee_[3*i+1];
            double dz = c[3*i+2] - coordinates_ref_wrapped_ee_[3*i+2];
            apply_min_image(dx, dy, dz);
            if (const double dr2 = dx*dx + dy*dy + dz*dz; dr2 > max_dr2) max_dr2 = dr2;
            if (max_dr2 > rebuild_disp2_ee_) return true; // early exit
        }
        return false;
    }

    void System::lj_switch_factors(double r, double& S, double& dSdr) const {
        if (!lj_use_switch_ || r <= lj_switch_) { S = 1.0; dSdr = 0.0; return; }
        if (r >= lj_cutoff_)                    { S = 0.0; dSdr = 0.0; return; }

        const double inv = 1.0 / (lj_cutoff_ - lj_switch_);
        const double x   = (r - lj_switch_) * inv;   // x in (0,1)

        const double x2 = x*x, x3 = x2*x, x4 = x2*x2, x5 = x4*x;

        // S = 1 - 10x^3 + 15x^4 - 6x^5  (same as OpenMM, just reordered)
        S = 1.0 - 10.0*x3 + 15.0*x4 - 6.0*x5;

        // dS/dx = -30x^2 + 60x^3 - 30x^4
        const double dSdx = (-30.0*x2 + 60.0*x3 - 30.0*x4);
        dSdr = dSdx * inv;
    }

    void System::build_nonbonded_cache()
    {
        type_.resize(natoms_);
        q_.resize(natoms_);

        // 1) type[] and q[]
        int maxType = -1;
        const auto& atoms = topology_.get_atom_list();
        for (size_t i = 0; i < natoms_; ++i) {
            const int ti = static_cast<int>(atoms[i].get_atom_type_index()) - 1; // 0-based
            type_[i] = ti;
            q_[i] = atoms[i].get_partial_charge();
            if (ti > maxType) maxType = ti;
        }

        nTypes_ = maxType + 1;

        // 2) flatten nbmatrix
        const auto& nb = topology_.get_nb_matrix(); // vector<vector<unsigned long>>
        nb_flat_.assign(static_cast<size_t>(nTypes_) * static_cast<size_t>(nTypes_), 0UL);

        for (int i = 0; i < nTypes_; ++i) {
            for (int j = 0; j < nTypes_; ++j) {
                nb_flat_[static_cast<size_t>(i) * nTypes_ + j] = nb[i][j];
            }
        }
    }

    void System::precompute_lj_energy_force_shift() {
        const auto& A   = topology_.get_lennard_jones_Acoefs_();
        const auto& B   = topology_.get_lennard_jones_Bcoefs_();
        const auto& A14 = topology_.get_lennard_jones_14_Acoefs_();
        const auto& B14 = topology_.get_lennard_jones_14_Bcoefs_();

        lj_Ucut_.resize(A.size());
        for (size_t p = 0; p < A.size(); ++p)
            lj_Ucut_[p] = NonBonded::LennardJones::CalculateEnergy(lj_cutoff2_, A[p], B[p]);

        lj_Ucut14_.resize(A14.size());
        for (size_t p = 0; p < A14.size(); ++p)
            lj_Ucut14_[p] = NonBonded::LennardJones::CalculateEnergy(lj_cutoff2_, A14[p], B14[p]);

        lj_Gcut_.resize(A.size());
        for (size_t p = 0; p < A.size(); ++p)
            lj_Gcut_[p] = NonBonded::LennardJones::CalculateGradient(lj_cutoff2_, A[p], B[p]);

        lj_Gcut14_.resize(A14.size());
        for (size_t p = 0; p < A14.size(); ++p)
            lj_Gcut14_[p] = NonBonded::LennardJones::CalculateGradient(lj_cutoff2_, A14[p], B14[p]);
    }

    void System::build_lj_pairlist()
    {
        if (!pbc_enabled_) return;
        lj_pairs_.clear();
        lj_pairs_.reserve(natoms_ * 60);

        const auto& coordinates = topology_.get_coordinates();
        head_.assign(nCells_lj_, -1);
        next_.assign(natoms_, -1);

        for (int i = 0; i < natoms_; ++i) // Linked list based cell-division
        {
            const double x = wrap(coordinates[3*i], boxLx_);
            const double y = wrap(coordinates[3*i + 1], boxLy_);
            const double z = wrap(coordinates[3*i + 2], boxLz_);

            const int ix = cell_number_along_dim(x, boxLx_, nCells_x_lj_);
            const int iy = cell_number_along_dim(y, boxLy_, nCells_y_lj_);
            const int iz = cell_number_along_dim(z, boxLz_, nCells_z_lj_);

            const int cell_ID = get_Cell_ID_lj(ix, iy, iz);
            next_[i] = head_[cell_ID];
            head_[cell_ID] = i;
        }

        for (int iz = 0; iz < nCells_z_lj_; ++iz)
        {
            for (int iy = 0; iy < nCells_y_lj_; ++iy)
            {
                for (int ix = 0; ix < nCells_x_lj_; ++ix)
                {
                    const int cell_ID_1 = get_Cell_ID_lj(ix, iy, iz);
                    std::array<int, 27> neigh{};
                    int m = 0;

                    for (int dz = -1; dz <= 1; ++dz)
                    {
                        for (int dy = -1; dy <= 1; ++dy)
                        {
                            for (int dx = -1; dx <= 1; ++dx)
                            {
                                const int jx = (ix + dx + nCells_x_lj_) % nCells_x_lj_;
                                const int jy = (iy + dy + nCells_y_lj_) % nCells_y_lj_;
                                const int jz = (iz + dz + nCells_z_lj_) % nCells_z_lj_;
                                neigh[m++] = get_Cell_ID_lj(jx, jy, jz);
                            }
                        }
                    }

                    std::sort(neigh.begin(), neigh.begin() + m);
                    m = static_cast<int>(std::unique(neigh.begin(), neigh.begin() + m) - neigh.begin());

                    for (int t = 0; t < m; ++t)
                    {
                        const int cell_ID_2 = neigh[t];
                        if (cell_ID_2 < cell_ID_1) continue;

                        for (int i = head_[cell_ID_1]; i != -1; i = next_[i])
                        {
                            for (int j = head_[cell_ID_2]; j != -1; j = next_[j]) {

                                if (cell_ID_1 == cell_ID_2 && j <= i) continue;

                                int a = i, b = j;
                                if (b < a) std::swap(a, b);

                                const auto& excl = topology_.get_atom_list()[a].get_excluded_atoms();

                                bool is14 = false;
                                if (std::ranges::binary_search(excl, b))
                                {
                                    is14 = topology_.is_14_pair(a, b);
                                    if (!is14) continue;
                                }

                                double dxv = coordinates[3*a]   - coordinates[3*b];
                                double dyv = coordinates[3*a+1] - coordinates[3*b+1];
                                double dzv = coordinates[3*a+2] - coordinates[3*b+2];
                                apply_min_image(dxv, dyv, dzv);

                                if (const double r2 = dxv*dxv + dyv*dyv + dzv*dzv; r2 <= lj_list_cutoff2_ && r2 > 1e-12)
                                {
                                    const int ti = type_[a];
                                    const int tj = type_[b];
                                    const unsigned long nb = nb_flat_[static_cast<size_t>(ti) * nTypes_ + tj];

                                    if (nb == 0) continue; // in the PBC while-loop version const uint32_t

                                    const auto p = static_cast<uint32_t>(nb - 1);
                                    lj_pairs_.push_back(NBPair{a, b, p, static_cast<uint8_t>(is14 ? 1 : 0)});
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void System::build_ee_pairlist()
    {
        if (!pbc_enabled_) return;
        ee_pairs_buf_.clear();
        ee_pairs_buf_.reserve(natoms_ * 60);

        const auto& coordinates = topology_.get_coordinates();
        head_.assign(nCells_ee_, -1);
        next_.assign(natoms_, -1);

        for (int i = 0; i < natoms_; ++i) // Linked list based cell-division
        {
            const double x = wrap(coordinates[3*i], boxLx_);
            const double y = wrap(coordinates[3*i + 1], boxLy_);
            const double z = wrap(coordinates[3*i + 2], boxLz_);

            const int ix = cell_number_along_dim(x, boxLx_, nCells_x_ee_);
            const int iy = cell_number_along_dim(y, boxLy_, nCells_y_ee_);
            const int iz = cell_number_along_dim(z, boxLz_, nCells_z_ee_);

            const int cell_ID = get_Cell_ID_ee(ix, iy, iz);
            next_[i] = head_[cell_ID];
            head_[cell_ID] = i;
        }

        for (int iz = 0; iz < nCells_z_ee_; ++iz)
        {
            for (int iy = 0; iy < nCells_y_ee_; ++iy)
            {
                for (int ix = 0; ix < nCells_x_ee_; ++ix)
                {
                    const int cell_ID_1 = get_Cell_ID_ee(ix, iy, iz);
                    std::array<int, 27> neigh{};
                    int m = 0;

                    for (int dz = -1; dz <= 1; ++dz)
                    {
                        for (int dy = -1; dy <= 1; ++dy)
                        {
                            for (int dx = -1; dx <= 1; ++dx)
                            {
                                const int jx = (ix + dx + nCells_x_ee_) % nCells_x_ee_;
                                const int jy = (iy + dy + nCells_y_ee_) % nCells_y_ee_;
                                const int jz = (iz + dz + nCells_z_ee_) % nCells_z_ee_;
                                neigh[m++] = get_Cell_ID_ee(jx, jy, jz);
                            }
                        }
                    }

                    std::sort(neigh.begin(), neigh.begin() + m);
                    m = static_cast<int>(std::unique(neigh.begin(), neigh.begin() + m) - neigh.begin());

                    for (int t = 0; t < m; ++t)
                    {
                        const int cell_ID_2 = neigh[t];
                        if (cell_ID_2 < cell_ID_1) continue;

                        for (int i = head_[cell_ID_1]; i != -1; i = next_[i])
                        {
                            for (int j = head_[cell_ID_2]; j != -1; j = next_[j]) {

                                if (cell_ID_1 == cell_ID_2 && j <= i) continue;

                                int a = i, b = j;
                                if (b < a) std::swap(a, b);

                                const auto& excl = topology_.get_atom_list()[a].get_excluded_atoms();

                                bool is14 = false;
                                if (std::ranges::binary_search(excl, b))
                                {
                                    is14 = topology_.is_14_pair(a, b);
                                    if (!is14) continue;
                                }

                                double dxv = coordinates[3*a]   - coordinates[3*b];
                                double dyv = coordinates[3*a+1] - coordinates[3*b+1];
                                double dzv = coordinates[3*a+2] - coordinates[3*b+2];
                                apply_min_image(dxv, dyv, dzv);

                                if (const double r2 = dxv*dxv + dyv*dyv + dzv*dzv; r2 <= ee_list_cutoff2_ && r2 > 1e-12)
                                {
                                    const int ti = type_[a];
                                    const int tj = type_[b];
                                    const unsigned long nb = nb_flat_[static_cast<size_t>(ti) * nTypes_ + tj];

                                    if (nb == 0) continue; // in the PBC while-loop version const uint32_t

                                    const auto p = static_cast<uint32_t>(nb - 1);
                                    ee_pairs_buf_.push_back(NBPair{a, b, p, static_cast<uint8_t>(is14 ? 1 : 0)});
                                }
                            }
                        }
                    }
                }
            }
        }
    }


}