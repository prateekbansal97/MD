#include "atom.h"
#include "CosineDihedral.h"
#include <cmath>

double CosineDihedral::return_energy(double angle)
{
    double energy = Dihedral_force_constant * (1.0 + std::cos(Dihedral_Periodicity * angle - Dihedral_Phase));
    return energy;
}