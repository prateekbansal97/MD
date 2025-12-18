//
// Created by Prateek Bansal on 12/9/25.
//

#ifndef MD_SYSTEM_H
#define MD_SYSTEM_H

#include <utility>
#include <vector>
#include <algorithm>

#include "../../include/AmberTopology/topology.h"

class System
{
public:
    explicit System(Topology  topology): topology(std::move(topology)) {init_forces();}

    Topology& get_topology() { return topology;}


    void init();
    void init_forces() {
        forces.assign(3 * topology.get_nAtoms(), 0.0);
    }

    // Zero out forces before every L-BFGS step
    void clear_forces() {
        std::fill(forces.begin(), forces.end(), 0.0);
    }

    void set_lj_cutoff(const double rc) { lj_cutoff_ = rc; lj_cutoff2_ = rc*rc; }
    void set_box(const double Lx, const double Ly, const double Lz) {
        boxLx_ = Lx; boxLy_ = Ly; boxLz_ = Lz;
        pbc_enabled_ = (Lx > 0 && Ly > 0 && Lz > 0);
    }

    void apply_min_image(double& dx, double& dy, double& dz) const;
    static double min_image_1d(double d, double L);
    void precompute_lj_energy_force_shift();



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

    bool pbc_enabled_ = false;
    double boxLx_ = 0.0, boxLy_ = 0.0, boxLz_ = 0.0;

    std::vector<double> lj_Ucut_;    // same indexing as A/B (nb-1)
    std::vector<double> lj_Ucut14_;
    std::vector<double> lj_Gcut_;    // same indexing as A/B (p = nb-1)
    std::vector<double> lj_Gcut14_;  // same indexing as A14/B14

};

#endif //MD_SYSTEM_H
