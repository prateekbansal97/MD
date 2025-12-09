
#include "../../include/AmberTopology/atom.h"
#include "../../include/AmberTopology/harmonicUB.h"

double HarmonicUB::return_energy(double distance)
{
    double force = UB_force_constant*(distance - UB_force_equil)*(distance - UB_force_equil);
    return force;
}

