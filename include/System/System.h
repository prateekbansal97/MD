//
// Created by Prateek Bansal on 12/9/25.
//

#ifndef MD_SYSTEM_H
#define MD_SYSTEM_H

#include <utility>

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
};

#endif //MD_SYSTEM_H
