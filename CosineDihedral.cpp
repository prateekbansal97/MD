#include "atom.h"
#include "CosineDihedral.h"

double CosineDihedral::return_energy(double angle)
{
    double energy = Dihedral_force_constant*(angle - Dihedral_Phase)*(angle - Dihedral_Phase);
    return energy;
}