//
// Created by Prateek Bansal on 12/13/25.
//

//18.2206899283247L

#include "../../include/AmberTopology/CoulombicEE.h"

double CoulombicEE::CalculateEnergy(double r, double chargeA, double chargeB, double epsilon)
{
    return ((chargeA * chargeB) / (r*epsilon)) * 18.2206899283247L * 18.2206899283247L;
}

double CoulombicEE::CalculateGradient(double r, double chargeA, double chargeB, double epsilon)
{
    double r3 = r * r * r;
    return -((chargeA * chargeB) / (r3*epsilon)) * 18.2206899283247L * 18.2206899283247L;
}