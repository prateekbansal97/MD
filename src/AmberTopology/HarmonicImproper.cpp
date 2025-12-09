
#include "../../include/AmberTopology/atom.h"
#include "../../include/AmberTopology/HarmonicImproper.h"
#include <cmath>

double HarmonicImproper::return_energy(double angle_radians)
{
    // CHARMM improper harmonic term: V = K_psi * (psi - psi0)^2
    // Returning the force-like value requested: 2 * K_psi * (psi - psi0)^2
    const double delta = angle_radians - IMP_Phase;
    return IMP_force_constant * delta * delta;
}
