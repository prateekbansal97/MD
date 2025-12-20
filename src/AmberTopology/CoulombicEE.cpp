//
// Created by Prateek Bansal on 12/13/25.
//

#include "AmberTopology/CoulombicEE.h"

double CoulombicEE::CalculateEnergy(const double r, const double chargeA, const double chargeB, const double epsilon)
{
    return ((chargeA * chargeB) / (r*epsilon)) * 18.2206899283247L * 18.2206899283247L;
}

double CoulombicEE::CalculateGradient(const double r, const double chargeA, const double chargeB, const double epsilon)
{
    const double r3 = r * r * r;
    return -((chargeA * chargeB) / (r3*epsilon)) * 18.2206899283247L * 18.2206899283247L;
}