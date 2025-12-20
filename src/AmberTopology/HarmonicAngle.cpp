
#include "AmberTopology/HarmonicAngle.h"

double HarmonicAngle::calculate_energy(const double angle) const
{
    const double energy = Angle_force_constant_*(angle - Angle_Equil_)*(angle - Angle_Equil_);
    return energy;
}
