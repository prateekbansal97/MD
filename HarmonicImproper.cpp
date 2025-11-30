
#include "atom.h"
#include "HarmonicImproper.h"

float HarmonicImproper::return_force(float distance)
{
    float force = IMP_force_constant*(distance - IMP_force_equil)*(distance - IMP_force_equil);
    return force;
}

