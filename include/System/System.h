//
// Created by Prateek Bansal on 12/9/25.
//

#ifndef MD_SYSTEM_H
#define MD_SYSTEM_H

#include <utility>
#include <vector>
#include <algorithm>
#include <cstdint>
// #include <array>
#include <ranges>

#include "../../include/AmberTopology/topology.h"


struct NBPair {
    int i;
    int j;
    uint32_t p;     // nb - 1
    uint8_t is14;   // 0 or 1
};

class System
{
public:
    explicit System(Topology  topology):
    topology(std::move(topology)),
    lj_cutoff2_(0.0), lj_skin_(0.0), lj_list_cutoff(0.0),
    lj_list_cutoff2(0.0), nCells_lj_(0), nCells_x_lj(0),
    nCells_y_lj(0), nCells_z_lj(0), lcell_x_lj(0.0),
    lcell_y_lj(0.0), lcell_z_lj(0.0),
    ee_cutoff2_(0.0), ee_list_cutoff(0.0),
    ee_list_cutoff2(0.0), nCells_ee_(0), nCells_x_ee(0),
    nCells_y_ee(0), nCells_z_ee(0), lcell_x_ee(0.0),
    lcell_y_ee(0.0), lcell_z_ee(0.0),
    natoms_(topology.get_nAtoms()) {init_forces();}

    Topology& get_topology() { return topology;}


    void init();
    void init_forces() {
        forces.assign(3 * topology.get_nAtoms(), 0.0);
    }

    // Zero out forces before every L-BFGS step
    void clear_forces() {
        std::ranges::fill(forces, 0.0);
    }

    void set_lj_cutoff(const double rc) { lj_cutoff_ = rc; lj_cutoff2_ = rc*rc; }
    void set_box(const double Lx, const double Ly, const double Lz) {
        boxLx_ = Lx; boxLy_ = Ly; boxLz_ = Lz;
        pbc_enabled_ = (Lx > 0 && Ly > 0 && Lz > 0);
        // if (pbc_enabled_) std::cout << "\nPBC Enabled! \n" ;
    }

    void set_lj_switch_radius(const double rswitch) {
        lj_switch_ = rswitch;
        lj_use_switch_ = (rswitch > 0.0);
        // optional safety:
        if (lj_use_switch_ && !(lj_switch_ < lj_cutoff_)) {
            throw std::runtime_error("lj_switch must be < lj_cutoff");
        }
    }

    void set_use_lj_switch(const bool use_lj) { lj_use_switch_ = use_lj; }

    void set_lj_skin_(const double skin)
    {
        lj_skin_ = skin;
        lj_list_cutoff = lj_skin_ + lj_cutoff_;
        lj_list_cutoff2 = lj_list_cutoff*lj_list_cutoff;
        // std::cout << "\nlj_list_cutoff2: " << lj_list_cutoff2 << std::endl;
    }

    void set_lj_list_cutoff(const double cutoff)
    {
        lj_list_cutoff = cutoff;
        lj_list_cutoff2 = cutoff*cutoff;
    }

    void set_ee_cutoff(const double cutoff)
    {
        ee_cutoff_ = cutoff;
        ee_cutoff2_ = cutoff*cutoff;
    }

    void set_ee_skin(const double eeskin)
    {
        ee_skin_ = eeskin;
        ee_list_cutoff = ee_skin_ + ee_cutoff_;
        ee_list_cutoff2 = ee_list_cutoff*ee_list_cutoff;
    }


    [[nodiscard]] double get_lj_list_cutoff() const {return lj_list_cutoff;}
    [[nodiscard]] double get_ee_list_cutoff() const {return ee_list_cutoff;}

    [[nodiscard]] double get_lj_skin() const {return lj_skin_;}

    void apply_min_image(double& dx, double& dy, double& dz) const;
    static double min_image_1d(double d, double L);
    void precompute_lj_energy_force_shift();
    static inline double wrap(double x, double Lx) ;
    [[nodiscard]] int get_Cell_ID_lj(int ix, int iy, int iz) const;
    [[nodiscard]] int get_Cell_ID_ee(int ix, int iy, int iz) const;

    [[nodiscard]] static int cell_number_along_dim(double x, double L, int ncelldim) ;
    void init_box();



    void calculate_energies();
    void calculate_bond_energy();
    void calculate_angle_energy();
    void calculate_dihedral_energy();
    void calculate_improper_energy();
    void calculate_UB_energy();
    void calculate_CMAP_energy();
    void calculate_LJ_energy();
    void calculate_LJ_energy_pairlist();
    void calculate_EE_energy();
    void calculate_EE_energy_pairlist();


    // This is the function your L-BFGS solver will call repeatedly
    void calculate_forces();
    void calculate_forces_bonds();
    void calculate_forces_UBbonds();
    void calculate_forces_angles();
    void calculate_forces_cosinedihedrals();
    void calculate_forces_harmonicImpropers();
    void calculate_forces_cmap();
    void calculate_forces_LJ();
    void calculate_forces_LJ_pairlist();
    void calculate_forces_EE();

    void build_nonbonded_cache();
    void build_lj_pairlist();
    void build_ee_pairlist();


    void setup_neighbor_rebuild_threshold_lj();
    void setup_neighbor_rebuild_threshold_ee();
    void copy_ref_coords_lj();
    void copy_ref_coords_ee();
    void setup_pairlists();
    void set_lj_switch(double rswitch);
    inline void lj_switch_factors(double r, double& S, double& dSdr) const;


    [[nodiscard]] bool ljpairs_needs_rebuild();
    [[nodiscard]] bool eepairs_needs_rebuild();



private:
    Topology topology;
    std::vector<double> forces;
    double bond_energies = 0.0;
    double angle_energy = 0.0;
    double dihedral_energy = 0.0;
    double urey_bradley_energy = 0.0;
    double improper_energy = 0.0;
    double CMAP_energy = 0.0;
    double LJ_energy = 0.0;
    double LJ_energy_pairlist = 0.0;
    double EE_energy = 0.0;
    double EE_energy_pairlist = 0.0;

    std::vector<int> type_;      // size N, 0-based type index
    std::vector<double> q_;      // size N, charges
    std::vector<unsigned long> nb_flat_;

    int nTypes_ = 0;

    double lj_cutoff_  = 10.0;                 // Angstrom
    double lj_cutoff2_ = lj_cutoff_ * lj_cutoff_;

    double lj_skin_;
    double lj_list_cutoff;
    double lj_list_cutoff2;

    double ee_cutoff_ = 10.0;
    double ee_cutoff2_ = ee_cutoff_ * ee_cutoff_;

    double ee_skin_ = 0.0;
    double ee_list_cutoff;
    double ee_list_cutoff2;


    bool pbc_enabled_ = false;
    double boxLx_ = 0.0, boxLy_ = 0.0, boxLz_ = 0.0;

    std::vector<double> lj_Ucut_;    // same indexing as A/B (nb-1)
    std::vector<double> lj_Ucut14_;
    std::vector<double> lj_Gcut_;    // same indexing as A/B (p = nb-1)
    std::vector<double> lj_Gcut14_;  // same indexing as A14/B14

    int nCells_lj_;
    int nCells_x_lj;
    int nCells_y_lj;
    int nCells_z_lj;

    double lcell_x_lj;
    double lcell_y_lj;
    double lcell_z_lj;

    int nCells_ee_;
    int nCells_x_ee;
    int nCells_y_ee;
    int nCells_z_ee;

    double lcell_x_ee;
    double lcell_y_ee;
    double lcell_z_ee;

    const unsigned long int natoms_;

    std::vector<NBPair> lj_pairs_;
    // std::vector<NBPair> ee_pairs_;

    // std::vector<NBPair> lj_pairs_;      // always built by build_lj_pairlist()
    std::vector<NBPair> ee_pairs_buf_;  // only built when thresholds differ

    std::span<const NBPair> ee_pairs_;

    std::vector<int> head;
    std::vector<int> next;

    std::vector<double> coordinates_ref_wrapped_lj;
    std::vector<double> coordinates_ref_wrapped_ee;
    double rebuild_disp2_lj = 0.0;
    double rebuild_disp2_ee = 0.0;

    double lj_switch_ = 0.0;     // rswitch
    bool   lj_use_switch_ = false;


};

#endif //MD_SYSTEM_H
