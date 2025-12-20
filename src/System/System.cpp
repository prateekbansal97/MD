//
// Created by Prateek Bansal on 12/9/25.
//

#include "System/System.h"
#include "System/Metrics.h"
#include "AmberTopology/Topology.h"
#include "AmberTopology/LennardJones.h"
#include "AmberTopology/Atom.h"
#include "AmberTopology/CoulombicEE.h"

#include <vector>
#include <cmath>
#include <algorithm>
#include <array>

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


void System::set_lj_switch(const double rswitch) {
    lj_switch_ = rswitch;
    lj_use_switch_ = (rswitch > 0.0);
    // optional safety:
    if (lj_use_switch_ && !(lj_switch_ < lj_cutoff_)) {
        throw std::runtime_error("lj_switch must be < lj_cutoff");
    }
}

inline void System::lj_switch_factors(const double r, double& S, double& dSdr) const {
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

inline double System::min_image_1d(double d, const double L) {
    // move into [-L/2, L/2]
    d -= L * std::round(d / L);
    return d;
}

inline void System::apply_min_image(double& dx, double& dy, double& dz) const {
    if (!pbc_enabled_) return;
    dx = min_image_1d(dx, boxLx_);
    dy = min_image_1d(dy, boxLy_);
    dz = min_image_1d(dz, boxLz_);
}

inline double System::wrap(double x, const double Lx)
{
    x -= Lx * std::floor(x/Lx);
    return x;
}

inline int System::get_Cell_ID_lj(const int ix, const int iy, const int iz) const
{
    return ix + nCells_x_lj_*(iy + nCells_y_lj_*iz);
}

inline int System::get_Cell_ID_ee(const int ix, const int iy, const int iz) const
{
    return ix + nCells_x_ee_*(iy + nCells_y_ee_*iz);
}

int System::cell_number_along_dim(const double x, const double L, const int nCelldim)
{
    // x is assumed wrapped into [0, L)
    int ix = static_cast<int>(x * nCelldim / L);   // 0..nCells-1
    if (ix >= nCelldim) ix = nCelldim - 1;           // safety for x ~ L due to fp
    if (ix < 0) ix = 0;
    return ix;
}

void System::init_box()
{
    if (!pbc_enabled_) return;
    nCells_x_lj_ = std::max(1, static_cast<int>(std::floor(boxLx_ / get_lj_list_cutoff())));
    nCells_y_lj_ = std::max(1, static_cast<int>(std::floor(boxLy_ / get_lj_list_cutoff())));
    nCells_z_lj_ = std::max(1, static_cast<int>(std::floor(boxLz_ / get_lj_list_cutoff())));

    lcell_x_lj_ = boxLx_ / nCells_x_lj_;
    lcell_y_lj_ = boxLy_ / nCells_y_lj_;
    lcell_z_lj_ = boxLz_ / nCells_z_lj_;

    nCells_lj_ = nCells_x_lj_ * nCells_y_lj_ * nCells_z_lj_;

    nCells_x_ee_ = std::max(1, static_cast<int>(std::floor(boxLx_ / get_ee_list_cutoff())));
    nCells_y_ee_ = std::max(1, static_cast<int>(std::floor(boxLy_ / get_ee_list_cutoff())));
    nCells_z_ee_ = std::max(1, static_cast<int>(std::floor(boxLz_ / get_ee_list_cutoff())));

    lcell_x_ee_ = boxLx_ / nCells_x_ee_;
    lcell_y_ee_ = boxLy_ / nCells_y_ee_;
    lcell_z_ee_ = boxLz_ / nCells_z_ee_;

    nCells_ee_ = nCells_x_ee_ * nCells_y_ee_ * nCells_z_ee_;
}

void System::build_nonbonded_cache()
{
    // const size_t N = topology.get_num_atoms();
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
        lj_Ucut_[p] = LennardJones::CalculateEnergy(lj_cutoff2_, A[p], B[p]);

    lj_Ucut14_.resize(A14.size());
    for (size_t p = 0; p < A14.size(); ++p)
        lj_Ucut14_[p] = LennardJones::CalculateEnergy(lj_cutoff2_, A14[p], B14[p]);

    lj_Gcut_.resize(A.size());
    for (size_t p = 0; p < A.size(); ++p)
        lj_Gcut_[p] = LennardJones::CalculateGradient(lj_cutoff2_, A[p], B[p]);

    lj_Gcut14_.resize(A14.size());
    for (size_t p = 0; p < A14.size(); ++p)
        lj_Gcut14_[p] = LennardJones::CalculateGradient(lj_cutoff2_, A14[p], B14[p]);
}

void System::build_lj_pairlist()
{
    if (!pbc_enabled_) return;
    // const size_t N = topology.get_num_atoms();
    lj_pairs_.clear();
    lj_pairs_.reserve(natoms_ * 60);

    const auto& coordinates = topology_.get_coordinates();
    // const double rlist = get_lj_list_cutoff();
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
    // const size_t N = topology.get_num_atoms();
    ee_pairs_buf_.clear();
    ee_pairs_buf_.reserve(natoms_ * 60);

    const auto& coordinates = topology_.get_coordinates();
    // const double rlist = get_ee_list_cutoff();
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
}

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

void System::calculate_LJ_energy()
{

    const std::vector<double>& coordinates = topology_.get_coordinates();
    std::vector<Atom>& atom_list = topology_.get_atom_list();

    // const size_t N = topology.get_num_atoms();

    const auto& A    = topology_.get_lennard_jones_Acoefs_();
    const auto& B    = topology_.get_lennard_jones_Bcoefs_();
    const auto& A14  = topology_.get_lennard_jones_14_Acoefs_();
    const auto& B14  = topology_.get_lennard_jones_14_Bcoefs_();

    LJ_energy_ = 0;

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

            // double distance_AB = distance(x1, y1, z1, x2, y2, z2);
            double dx = x1 - x2;
            double dy = y1 - y2;
            double dz = z1 - z2;
            apply_min_image(dx,dy,dz);
            const double r2 = dx*dx + dy*dy + dz*dz;

            if (r2 > lj_cutoff2_) continue;
            if (r2 < 1e-12) continue;

            const int ti = type_[atomAIndex];
            const int tj = type_[atomBIndex];
            const unsigned long nb = nb_flat_[static_cast<size_t>(ti) * nTypes_ + tj];

            if (nb == 0) continue;

            double Aij = 0, Bij = 0;
            // const auto p = nb - 1;
            const auto p = static_cast<size_t>(nb - 1);

            const double r = std::sqrt(r2);

            if (!is14) {
                Aij = A[p];
                Bij = B[p];
            } else {
                Aij = A14[p];
                Bij = B14[p];
            }

            // const double energy = LennardJones::CalculateEnergy(r2, Aij, Bij);
            double energy = LennardJones::CalculateEnergy(r2, Aij, Bij);
            energy -= (!is14 ? lj_Ucut_[p] : lj_Ucut14_[p]);
            const double gcut = (!is14 ? lj_Gcut_[p] : lj_Gcut14_[p]);
            energy += (r - lj_cutoff_) * (gcut * lj_cutoff_);
            LJ_energy_ += energy;
        }
    }
}

void System::calculate_LJ_energy_pairlist()
{

    LJ_energy_pairlist_ = 0;
    const std::vector<double>& coordinates = topology_.get_coordinates();
    // std::vector<Atom>& atom_list = topology.get_atoms();

    // const size_t N = topology.get_num_atoms();

    const auto& A    = topology_.get_lennard_jones_Acoefs_();
    const auto& B    = topology_.get_lennard_jones_Bcoefs_();
    const auto& A14  = topology_.get_lennard_jones_14_Acoefs_();
    const auto& B14  = topology_.get_lennard_jones_14_Bcoefs_();

    for (auto& [atomAIndex, atomBIndex, p, is14]: lj_pairs_)
    {

        const double x1 = coordinates[3*atomAIndex], y1 = coordinates[3*atomAIndex + 1], z1 = coordinates[3*atomAIndex + 2];
        const double x2 = coordinates[3*atomBIndex], y2 = coordinates[3*atomBIndex + 1], z2 = coordinates[3*atomBIndex + 2];

        // double distance_AB = distance(x1, y1, z1, x2, y2, z2);
        double dx = x1 - x2;
        double dy = y1 - y2;
        double dz = z1 - z2;
        apply_min_image(dx,dy,dz);
        const double r2 = dx*dx + dy*dy + dz*dz;

        if (r2 > lj_cutoff2_) continue;
        if (r2 < 1e-12) continue;

        double Aij = 0, Bij = 0;

        const double r = std::sqrt(r2);

        if (!is14) {
            Aij = A[p];
            Bij = B[p];
        } else {
            Aij = A14[p];
            Bij = B14[p];
        }

        // const double energy = LennardJones::CalculateEnergy(r2, Aij, Bij);
        double energy = LennardJones::CalculateEnergy(r2, Aij, Bij);

        if (lj_use_switch_) {
            double S, dSdr;
            lj_switch_factors(r, S, dSdr); // for r<=rswitch this should return S=1
            energy *= S;
        } else {
            energy -= (!is14 ? lj_Ucut_[p] : lj_Ucut14_[p]);
            const double gcut = (!is14 ? lj_Gcut_[p] : lj_Gcut14_[p]);
            energy += (r - lj_cutoff_) * (gcut * lj_cutoff_);
        }

        LJ_energy_pairlist_ += energy;
    }
}

void System::calculate_EE_energy()
{
    const std::vector<double>& coordinates = topology_.get_coordinates();
    std::vector<Atom>& atom_list = topology_.get_atom_list();
    // const size_t N = topology.get_num_atoms();

    EE_energy_ = 0;
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

            // double distance_AB = distance(x1, y1, z1, x2, y2, z2);
            double dx = x1 - x2;
            double dy = y1 - y2;
            double dz = z1 - z2;
            apply_min_image(dx,dy,dz);
            const double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 < 1e-12) continue;
            if (r2 > ee_cutoff2_) continue;

            const double r = std::sqrt(r2);

            double scale = 1.0;
            if (is14) scale = topology_.get_scee_scale_for_pair(static_cast<int>(atomAIndex), static_cast<int>(atomBIndex));

            const double energy = (1/scale) * CoulombicEE::CalculateEnergy(r, q_[atomAIndex], q_[atomBIndex], 1);

            EE_energy_ += energy;
        }
    }
}

void System::calculate_EE_energy_pairlist()
{
    const std::vector<double>& coordinates = topology_.get_coordinates();
    // std::vector<Atom>& atom_list = topology.get_atoms();
    // const size_t N = topology.get_num_atoms();

    EE_energy_pairlist_ = 0;
    for (auto& [atomAIndex, atomBIndex, p, is14]: ee_pairs_)
    {
        // Atom& atomA = atom_list[atomAIndex];
        // const std::vector<int>& excl = atomA.get_excluded_atoms();

        const double x1 = coordinates[3*atomAIndex], y1 = coordinates[3*atomAIndex + 1], z1 = coordinates[3*atomAIndex + 2];
        const double x2 = coordinates[3*atomBIndex], y2 = coordinates[3*atomBIndex + 1], z2 = coordinates[3*atomBIndex + 2];

        // double distance_AB = distance(x1, y1, z1, x2, y2, z2);
        double dx = x1 - x2;
        double dy = y1 - y2;
        double dz = z1 - z2;
        apply_min_image(dx,dy,dz);
        const double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 < 1e-12) continue;
        if (r2 > ee_cutoff2_) continue;

        const double r = std::sqrt(r2);

        double scale = 1.0;
        if (is14) scale = topology_.get_scee_scale_for_pair(static_cast<int>(atomAIndex), static_cast<int>(atomBIndex));

        const double energy = (1/scale) * CoulombicEE::CalculateEnergy(r, q_[atomAIndex], q_[atomBIndex], 1);

        EE_energy_pairlist_ += energy;
    }
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
    calculate_forces_EE();
}

void System::calculate_forces_bonds() {
    const std::vector<double>& coordinates = topology_.get_coordinates();

    for (auto& bond: topology_.get_harmonic_bonds())
    {
        const int atomAIndex = bond.get_atomA_index();
        const int atomBIndex = bond.get_atomB_index();
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

        forces_[3*atomAIndex] += fax; forces_[3*atomAIndex+1] += fay; forces_[3*atomAIndex+2] += faz;
        forces_[3*atomBIndex] -= fax; forces_[3*atomBIndex+1] -= fay; forces_[3*atomBIndex+2] -= faz;
    }
}

void System::calculate_forces_UBbonds() {
    // clear_forces(); // Reset gradients to 0
    const std::vector<double>& coordinates = topology_.get_coordinates();

    for (auto& bond: topology_.get_harmonic_UBs())
    {
        const int atomAIndex = bond.get_atomA_index();
        const int atomBIndex = bond.get_atomB_index();
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

        forces_[3*atomAIndex] += fax; forces_[3*atomAIndex+1] += fay; forces_[3*atomAIndex+2] += faz;
        forces_[3*atomBIndex] += -fax; forces_[3*atomBIndex+1] += -fay; forces_[3*atomBIndex+2] += -faz;
    }
}

void System::calculate_forces_angles() {
//    clear_forces(); // Reset gradients to 0
    const std::vector<double>& coordinates = topology_.get_coordinates();

    for (auto& ang: topology_.get_harmonic_angles())
    {
        const int atomAIndex = ang.get_atomA_index();
        const int atomBIndex = ang.get_atomB_index();
        const int atomCIndex = ang.get_atomC_index();

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
        const double theta = Metrics::angle(x1, y1, z1, x2, y2, z2, x3, y3, z3);

        const double fmag = -2 * k * (theta - theta0);
        const double fax = fmag/mod_ab*pa_x;
        const double fay = fmag/mod_ab*pa_y;
        const double faz = fmag/mod_ab*pa_z;

        const double fcx = fmag/mod_bc*pc_x;
        const double fcy = fmag/mod_bc*pc_y;
        const double fcz = fmag/mod_bc*pc_z;

        forces_[3*atomAIndex] += fax; forces_[3*atomAIndex+1] += fay; forces_[3*atomAIndex+2] += faz;
        forces_[3*atomCIndex] += fcx; forces_[3*atomCIndex+1] += fcy; forces_[3*atomCIndex+2] += fcz;
        forces_[3*atomBIndex] += -fax-fcx; forces_[3*atomBIndex+1] += -fay-fcy; forces_[3*atomBIndex+2] += -faz-fcz;

    }
}

void System::calculate_forces_cosinedihedrals() {
    const std::vector<double>& coordinates = topology_.get_coordinates();

    for (auto& dih: topology_.get_cosine_dihedrals()) {
        const int atomAIndex = dih.get_atomA_index();
        const int atomBIndex = dih.get_atomB_index();
        const int atomCIndex = dih.get_atomC_index();
        const int atomDIndex = dih.get_atomD_index();

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
        const double phi = Metrics::dihedral(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4);
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

        forces_[3*atomAIndex] += fax; forces_[3*atomAIndex+1] += fay; forces_[3*atomAIndex+2] += faz;
        forces_[3*atomBIndex] += fbx; forces_[3*atomBIndex+1] += fby; forces_[3*atomBIndex+2] += fbz;
        forces_[3*atomCIndex] += fcx; forces_[3*atomCIndex+1] += fcy; forces_[3*atomCIndex+2] += fcz;
        forces_[3*atomDIndex] += fdx; forces_[3*atomDIndex+1] += fdy; forces_[3*atomDIndex+2] += fdz;
    }
}

void System::calculate_forces_harmonicImpropers() {
    const std::vector<double>& coordinates = topology_.get_coordinates();

    for (auto& dih: topology_.get_harmonic_impropers()) {
        const int atomAIndex = dih.get_atomA_index();
        const int atomBIndex = dih.get_atomB_index();
        const int atomCIndex = dih.get_atomC_index();
        const int atomDIndex = dih.get_atomD_index();

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
        const double psi = Metrics::dihedral(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4);

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

        forces_[3*atomAIndex] += fax; forces_[3*atomAIndex+1] += fay; forces_[3*atomAIndex+2] += faz;
        forces_[3*atomBIndex] += fbx; forces_[3*atomBIndex+1] += fby; forces_[3*atomBIndex+2] += fbz;
        forces_[3*atomCIndex] += fcx; forces_[3*atomCIndex+1] += fcy; forces_[3*atomCIndex+2] += fcz;
        forces_[3*atomDIndex] += fdx; forces_[3*atomDIndex+1] += fdy; forces_[3*atomDIndex+2] += fdz;
    }
}

void System::calculate_forces_cmap()
{
    const std::vector<double>& coordinates = topology_.get_coordinates();
    for (auto& cmap : topology_.get_cmaps()) {

        const int atomAIndex = cmap.get_atomA_index();
        const int atomBIndex = cmap.get_atomB_index();
        const int atomCIndex = cmap.get_atomC_index();
        const int atomDIndex = cmap.get_atomD_index();
        const int atomEIndex = cmap.get_atomE_index();

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


        const double phi = Metrics::dihedral(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4);
        const double psi = Metrics::dihedral(x2, y2, z2, x3, y3, z3, x4, y4, z4, x5, y5, z5);

        // Calculate slopes (gradients)
        const int set_index = cmap.get_parameter_set();
        const std::vector<int>& resolutions = topology_.get_charmm_cmap_resolutions();
        const int resolution = resolutions[set_index - 1];
        const std::vector<double>& full_grid = topology_.get_cmap_grid_data(set_index);

        auto [fst, snd] = cmap.return_gradient_bicubic(phi, psi, resolution, topology_.get_Cmap_Coefficient_Matrix_bicubic_spline(), full_grid);
        const double dEdPhi = fst;
        const double dEdPsi = snd;



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


        forces_[3*atomAIndex] += fax_phi; forces_[3*atomAIndex+1] += fay_phi; forces_[3*atomAIndex+2] += faz_phi;
        forces_[3*atomBIndex] += fbx_phi; forces_[3*atomBIndex+1] += fby_phi; forces_[3*atomBIndex+2] += fbz_phi;
        forces_[3*atomCIndex] += fcx_phi; forces_[3*atomCIndex+1] += fcy_phi; forces_[3*atomCIndex+2] += fcz_phi;
        forces_[3*atomDIndex] += fdx_phi; forces_[3*atomDIndex+1] += fdy_phi; forces_[3*atomDIndex+2] += fdz_phi;

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

        forces_[3*atomBIndex] += fbx_psi; forces_[3*atomBIndex+1] += fby_psi; forces_[3*atomBIndex+2] += fbz_psi;
        forces_[3*atomCIndex] += fcx_psi; forces_[3*atomCIndex+1] += fcy_psi; forces_[3*atomCIndex+2] += fcz_psi;
        forces_[3*atomDIndex] += fdx_psi; forces_[3*atomDIndex+1] += fdy_psi; forces_[3*atomDIndex+2] += fdz_psi;
        forces_[3*atomEIndex] += fex_psi; forces_[3*atomEIndex+1] += fey_psi; forces_[3*atomEIndex+2] += fez_psi;
    }
}

void System::calculate_forces_LJ()
{
    const std::vector<double>& coordinates = topology_.get_coordinates();
    std::vector<Atom>& atom_list = topology_.get_atom_list();

    // const size_t N = topology.get_num_atoms();
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

            // double distance_AB = distance(x1, y1, z1, x2, y2, z2);
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
            // const auto p = nb - 1;
            const auto p = static_cast<size_t>(nb - 1);

            if (!is14) {
                Aij = A[p];
                Bij = B[p];
            } else {
                Aij = A14[p];
                Bij = B14[p];
            }

            // const double gradient = LennardJones::CalculateGradient(r2, Aij, Bij);

            double gradient = LennardJones::CalculateGradient(r2, Aij, Bij);
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

        // double distance_AB = distance(x1, y1, z1, x2, y2, z2);
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


        double gradient = LennardJones::CalculateGradient(r2, Aij, Bij);
        // const double gcut = (!is14 ? lj_Gcut_[p] : lj_Gcut14_[p]);
        // gradient -= gcut * (lj_cutoff_ / r);

        if (lj_use_switch_) {
            double S, dSdr;
            lj_switch_factors(r, S, dSdr);
            const double U = LennardJones::CalculateEnergy(r2, Aij, Bij);
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
    // const size_t N = topology.get_num_atoms();


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

            // double distance_AB = distance(x1, y1, z1, x2, y2, z2);
            double dx = x1 - x2;
            double dy = y1 - y2;
            double dz = z1 - z2;
            apply_min_image(dx,dy,dz);
            const double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 < 1e-12) continue;
            const double r = std::sqrt(r2);

            double scale = 1.0;
            if (is14) scale = topology_.get_scee_scale_for_pair(static_cast<int>(atomAIndex), static_cast<int>(atomBIndex));

            const double gradient = (1/scale) * CoulombicEE::CalculateGradient(r, q_[atomAIndex], q_[atomBIndex], 1);
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





