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


struct LJPair {
    int i;
    int j;
    uint32_t p;     // nb - 1
    uint8_t is14;   // 0 or 1
};

class System
{
public:
    explicit System(Topology  topology): topology(std::move(topology)), lj_cutoff2_(0.0), lj_skin_(0.0), lj_list_cutoff(0.0),
    lj_list_cutoff2(0.0), nCells_(0), nCells_x(0), nCells_y(0), nCells_z(0), lcell_x(0.0), lcell_y(0.0), lcell_z(0.0), natoms_(topology.get_nAtoms()) {init_forces();}

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
    }
    void set_lj_skin_(const double skin)
    {
        lj_skin_ = skin;
        lj_list_cutoff = lj_skin_ + lj_cutoff_;
        lj_list_cutoff2 = lj_list_cutoff*lj_list_cutoff;
    }

    void set_lj_list_cutoff(const double cutoff)
    {
        lj_list_cutoff = cutoff;
        lj_list_cutoff2 = cutoff*cutoff;
    }
    [[nodiscard]] double get_lj_list_cutoff() const {return lj_list_cutoff;}
    [[nodiscard]] double get_lj_skin() const {return lj_skin_;}

    void apply_min_image(double& dx, double& dy, double& dz) const;
    static double min_image_1d(double d, double L);
    void precompute_lj_energy_force_shift();
    static inline double wrap(double x, double Lx) ;
    [[nodiscard]] int get_Cell_ID(int ix, int iy, int iz) const;
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
    void calculate_EE_energy();


    // This is the function your L-BFGS solver will call repeatedly
    void calculate_forces();
    void calculate_forces_bonds();
    void calculate_forces_UBbonds();
    void calculate_forces_angles();
    void calculate_forces_cosinedihedrals();
    void calculate_forces_harmonicImpropers();
    void calculate_forces_cmap();
    void calculate_forces_LJ();
    void calculate_forces_EE();

    void build_nonbonded_cache();
    void build_lj_pairlist();


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
    double EE_energy = 0.0;

    std::vector<int> type_;      // size N, 0-based type index
    std::vector<double> q_;      // size N, charges
    std::vector<unsigned long> nb_flat_;

    int nTypes_ = 0;
    double lj_cutoff_  = 10.0;                 // Angstrom
    double lj_cutoff2_ = lj_cutoff_ * lj_cutoff_;

    double lj_skin_;
    double lj_list_cutoff;
    double lj_list_cutoff2;

    bool pbc_enabled_ = false;
    double boxLx_ = 0.0, boxLy_ = 0.0, boxLz_ = 0.0;

    std::vector<double> lj_Ucut_;    // same indexing as A/B (nb-1)
    std::vector<double> lj_Ucut14_;
    std::vector<double> lj_Gcut_;    // same indexing as A/B (p = nb-1)
    std::vector<double> lj_Gcut14_;  // same indexing as A14/B14

    int nCells_;
    int nCells_x;
    int nCells_y;
    int nCells_z;

    double lcell_x;
    double lcell_y;
    double lcell_z;

    const unsigned long int natoms_;
    std::vector<LJPair> lj_pairs_;
    std::vector<int> head;
    std::vector<int> next;

};

#endif //MD_SYSTEM_H
