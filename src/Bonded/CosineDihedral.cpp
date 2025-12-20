#include "Bonded/CosineDihedral.h"
#include <cmath>

namespace md::Bonded
{
    double CosineDihedral::calculate_energy(const double dihedra) const
    {
        const double d_energy = Dihedral_force_constant_ * (1.0 + std::cos(Dihedral_Periodicity_ * dihedra - Dihedral_Phase_));
        return d_energy;
    }
}