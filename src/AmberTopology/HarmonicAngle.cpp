
// #include "../../include/AmberTopology/atom.h"
#include "../../include/AmberTopology/HarmonicAngle.h"

double HarmonicAngle::return_energy(const double angle_) const
{
    const double energy_ = Angle_force_constant*(angle_ - Angle_Equil)*(angle_ - Angle_Equil);
    return energy_;
}
