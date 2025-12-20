
#include "AmberTopology/harmonicUB.h"

double HarmonicUB::return_energy(const double distance_) const
{
    const double force = UB_force_constant*(distance_ - UB_force_equil)*(distance_ - UB_force_equil);
    return force;
}

