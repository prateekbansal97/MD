
// #include "../../include/AmberTopology/atom.h"
#include "../../include/AmberTopology/HarmonicBond.h"

double HarmonicBond::return_energy(const double distance_) const
{
    const double force = Bond_force_constant*(distance_ - Bond_Equil)*(distance_ - Bond_Equil);
    return force;
}
