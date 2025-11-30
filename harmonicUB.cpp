
#include "atom.h"
#include "harmonicUB.h"

float HarmonicUB::return_force(float distance)
{
    float force = UB_force_constant*(distance - UB_force_equil)*(distance - UB_force_equil);
    return force;
}

