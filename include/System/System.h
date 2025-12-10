//
// Created by Prateek Bansal on 12/9/25.
//

#ifndef MD_SYSTEM_H
#define MD_SYSTEM_H

#include "../../include/AmberTopology/topology.h"

class System
{
public:
    System(Topology topology): topology(topology) {}
    Topology& get_topology() { return topology;}


    void init();
    double distance(double x1, double y1, double z1, double x2, double y2, double z2);
    double angle(double x1, double y1, double z1, double x2, double y2, double z2,  double x3, double y3, double z3);
    double dihedral(double x1, double y1, double z1, double x2, double y2, double z2,
                    double x3, double y3, double z3, double x4, double y4, double z4);


private:
    Topology topology;
};

#endif //MD_SYSTEM_H
