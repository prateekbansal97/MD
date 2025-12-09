
#include "../../include/AmberTopology/atom.h"
#include "../../include/AmberTopology/HarmonicAngle.h"

double HarmonicAngle::return_energy(double angle)
{
    double energy = Angle_force_constant*(angle - Angle_Equil)*(angle - Angle_Equil);
    return energy;
}