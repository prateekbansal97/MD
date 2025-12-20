
#include "Bonded/HarmonicBond.h"

namespace md::Bonded
{
    double HarmonicBond::calculate_energy(const double distance) const
    {
        const double force = Bond_force_constant_*(distance - Bond_Equil_)*(distance - Bond_Equil_);
        return force;
    }
}
