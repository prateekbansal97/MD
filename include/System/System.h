//
// Created by Prateek Bansal on 12/9/25.
//

#ifndef MD_SYSTEM_H
#define MD_SYSTEM_H

#include "../../include/AmberTopology/topology.h"

class System
{
public:
    System(Topology topology): topology(topology) {init_forces();}
    Topology& get_topology() { return topology;}


    void init();
    double distance(double x1, double y1, double z1, double x2, double y2, double z2);
    double angle(double x1, double y1, double z1, double x2, double y2, double z2,  double x3, double y3, double z3);
    double dihedral(double x1, double y1, double z1, double x2, double y2, double z2,
                    double x3, double y3, double z3, double x4, double y4, double z4);

    void init_forces() {
        forces.assign(3 * topology.get_nAtoms(), 0.0);
    }

    // Zero out forces before every L-BFGS step
    void clear_forces() {
        std::fill(forces.begin(), forces.end(), 0.0);
    }

    // This is the function your L-BFGS solver will call repeatedly
    void calculate_forces_bonds();
    void calculate_forces_UBbonds();
    void calculate_forces_angles();
    void calculate_forces_cosinedihedrals();

private:
    Topology topology;
    std::vector<double> forces;
};

#endif //MD_SYSTEM_H
