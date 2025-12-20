
#include "Bonded/harmonicUB.h"

namespace md::Bonded
{
    double HarmonicUB::calculate_energy(const double distance) const
    {
        const double force = UB_force_constant_*(distance - UB_force_equil_)*(distance - UB_force_equil_);
        return force;
    }
} //namespace md

