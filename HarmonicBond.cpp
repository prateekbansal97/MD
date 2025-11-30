
#include "atom.h"
#include "HarmonicBond.h"

double HarmonicBond::return_energy(double distance)
{
    double force = Bond_force_constant*(distance - Bond_Equil)*(distance - Bond_Equil);
    return force;
}