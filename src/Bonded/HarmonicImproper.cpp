
#include "Bonded/HarmonicImproper.h"

namespace md::Bonded
{
    double HarmonicImproper::calculate_energy(const double angle_radians) const
    {
        // CHARMM improper harmonic term: V = K_psi * (psi - psi0)^2
        // Returning the force-like value requested: 2 * K_psi * (psi - psi0)^2
        const double delta = angle_radians - IMP_Phase_;
        return IMP_force_constant_ * delta * delta;
    }
} // namespace md
