#include "AmberTopology/CosineDihedral.h"
#include <cmath>

double CosineDihedral::return_energy(const double dihedra) const
{
    const double d_energy = Dihedral_force_constant * (1.0 + std::cos(Dihedral_Periodicity * dihedra - Dihedral_Phase));
    return d_energy;
}
