//
// Created by Prateek Bansal on 12/9/25.
//

#ifndef MD_SYSTEM_H
#define MD_SYSTEM_H

#include <utility>
#include <vector>
#include <algorithm>
#include <span>

#include "AmberTopology/Topology.h"

namespace md
{
    struct NBPair {
        int i;
        int j;
        uint32_t p;     // nb - 1
        uint8_t is14;   // 0 or 1
    };

    class System
    {
    public:
        void init();

        
        explicit System(Topology  topology):
        topology_(std::move(topology)),
        lj_cutoff2_(0.0), lj_skin_(0.0), lj_list_cutoff_(0.0),
        lj_list_cutoff2_(0.0), ee_cutoff2_(0.0), ee_list_cutoff_(0.0),
        ee_list_cutoff2_(0.0), nCells_lj_(0), nCells_x_lj_(0),
        nCells_y_lj_(0), nCells_z_lj_(0),
        lcell_x_lj_(0.0), lcell_y_lj_(0.0),
        lcell_z_lj_(0.0), nCells_ee_(0), nCells_x_ee_(0),
        nCells_y_ee_(0), nCells_z_ee_(0), lcell_x_ee_(0.0),
        lcell_y_ee_(0.0), lcell_z_ee_(0.0),
        natoms_(topology.get_nAtoms()) {init_forces();}

        [[nodiscard]] const md::Topology& get_topology() const { return topology_;}

        void set_lj_cutoff(const double rc) { lj_cutoff_ = rc; lj_cutoff2_ = rc*rc; }
        void set_box(double Lx, double Ly, double Lz);
        void set_lj_switch_radius(const double rswitch);
        void set_use_lj_switch(const bool use_lj) { lj_use_switch_ = use_lj; }
        void set_lj_skin_(double skin);
        void set_lj_list_cutoff(double cutoff);
        void set_ee_cutoff(double cutoff);
        void set_ee_skin(double eeskin);

        [[nodiscard]] double get_lj_list_cutoff() const {return lj_list_cutoff_;}
        [[nodiscard]] double get_ee_list_cutoff() const {return ee_list_cutoff_;}
        [[nodiscard]] double get_lj_skin() const {return lj_skin_;}
        [[nodiscard]] int get_Cell_ID_lj(int ix, int iy, int iz) const;
        [[nodiscard]] int get_Cell_ID_ee(int ix, int iy, int iz) const;
        [[nodiscard]] static int cell_number_along_dim(double x, double L, int ncelldim);
        [[nodiscard]] bool ljpairs_needs_rebuild();
        [[nodiscard]] bool eepairs_needs_rebuild();



    private:

        void init_forces() {
            forces_.assign(3 * topology_.get_nAtoms(), 0.0);
        }

        // Zero out forces before every L-BFGS step
        void clear_forces() {
            std::ranges::fill(forces_, 0.0);
        }


        void apply_min_image(double& dx, double& dy, double& dz) const;
        static double min_image_1d(double d, double L);
        void precompute_lj_energy_force_shift();
        static double wrap(double x, double Lx) ;


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
        void lj_switch_factors(double r, double& S, double& dSdr) const;


        Topology topology_;
        std::vector<double> forces_;
        double bond_energies_ = 0.0;
        double angle_energy_ = 0.0;
        double dihedral_energy_ = 0.0;
        double urey_bradley_energy_ = 0.0;
        double improper_energy_ = 0.0;
        double CMAP_energy_ = 0.0;
        double LJ_energy_ = 0.0;
        double LJ_energy_pairlist_ = 0.0;
        double EE_energy_ = 0.0;
        double EE_energy_pairlist_ = 0.0;

        std::vector<int> type_;      // size N, 0-based type index
        std::vector<double> q_;      // size N, charges
        std::vector<unsigned long> nb_flat_;

        int nTypes_ = 0;

        double lj_cutoff_  = 10.0;                 // Angstrom
        double lj_cutoff2_ = lj_cutoff_ * lj_cutoff_;

        double lj_skin_;
        double lj_list_cutoff_;
        double lj_list_cutoff2_;

        double ee_cutoff_ = 10.0;
        double ee_cutoff2_ = ee_cutoff_ * ee_cutoff_;

        double ee_skin_ = 0.0;
        double ee_list_cutoff_;
        double ee_list_cutoff2_;


        bool pbc_enabled_ = false;
        double boxLx_ = 0.0, boxLy_ = 0.0, boxLz_ = 0.0;

        std::vector<double> lj_Ucut_;    // same indexing as A/B (nb-1)
        std::vector<double> lj_Ucut14_;
        std::vector<double> lj_Gcut_;    // same indexing as A/B (p = nb-1)
        std::vector<double> lj_Gcut14_;  // same indexing as A14/B14

        int nCells_lj_;
        int nCells_x_lj_;
        int nCells_y_lj_;
        int nCells_z_lj_;

        double lcell_x_lj_;
        double lcell_y_lj_;
        double lcell_z_lj_;

        int nCells_ee_;
        int nCells_x_ee_;
        int nCells_y_ee_;
        int nCells_z_ee_;

        double lcell_x_ee_;
        double lcell_y_ee_;
        double lcell_z_ee_;

        const unsigned long int natoms_;

        std::vector<NBPair> lj_pairs_;

        std::vector<NBPair> ee_pairs_buf_;  // only built when thresholds differ
        std::span<const NBPair> ee_pairs_;

        std::vector<int> head_;
        std::vector<int> next_;

        std::vector<double> coordinates_ref_wrapped_lj_;
        std::vector<double> coordinates_ref_wrapped_ee_;
        double rebuild_disp2_lj_ = 0.0;
        double rebuild_disp2_ee_ = 0.0;

        double lj_switch_ = 0.0;     // rswitch
        bool   lj_use_switch_ = false;

    };
} // namespace md

#endif //MD_SYSTEM_H
