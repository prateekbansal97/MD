//
// Created by Prateek Bansal on 12/13/25.
//

#include "NonBonded/CoulombicEE.h"
#include <cmath>

namespace md::NonBonded
{
    double CoulombicEE::CalculateEnergy(const double r, const double chargeA, const double chargeB, const double epsilon)
    {
        return ((chargeA * chargeB) / (r*epsilon)) * 18.2206899283247L * 18.2206899283247L;
    }

    double CoulombicEE::CalculateGradient(const double r, const double chargeA, const double chargeB, const double epsilon)
    {
        const double r3 = r * r * r;
        return ((chargeA * chargeB) / (r3*epsilon)) * 18.2206899283247L * 18.2206899283247L;
    }

    double CoulombicEE::CalculateEnergyCutoff(const double r, const double rcutoff, const double chargeA, const double chargeB, const double epsilon)
    {

        const double r2 = r * r;
        const double r_inv = 1.0 / r;
        const double rcutoff_inv = 1.0 / rcutoff;
        const double krf = rcutoff_inv * rcutoff_inv * rcutoff_inv * ((epsilon - 1)/(2*epsilon + 1));
        const double crf = rcutoff_inv * ((3*epsilon)/(2*epsilon + 1));
        const double cutoffenergy = (((chargeA * chargeB) / (epsilon)) * 18.2206899283247L * 18.2206899283247L)*(r_inv + krf*r2 - crf);
        return cutoffenergy;
    }

    double CoulombicEE::CalculateEnergyEwaldDirectTerm(const double r, const double chargeA, const double chargeB, const double ewald_alpha)
    {
        return ((chargeA * chargeB) / (r)) * std::erfc(ewald_alpha * r) * 18.2206899283247L * 18.2206899283247L;
    }

    double CoulombicEE::CalculateEnergyEwaldDirectTerm_14(const double r, const double chargeA, const double chargeB, const double ewald_alpha, const double scale)
    {

        const double energy_r = (chargeA * chargeB) * 18.2206899283247L * 18.2206899283247L;
        const double energy = (((1)/(scale*r)) - ((std::erf(r))/(r)))*energy_r;
        return energy;
    }

} // namespace md