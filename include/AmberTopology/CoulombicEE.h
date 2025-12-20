//
// Created by Prateek Bansal on 12/13/25.
//

#ifndef MD_COULOMBICEE_H
#define MD_COULOMBICEE_H

class CoulombicEE
{
private:
    static double CalculateEnergy(double r, double chargeA, double chargeB, double epsilon);
    static double CalculateGradient(double r, double chargeA, double chargeB, double epsilon);
    friend class System;
};

#endif //MD_COULOMBICEE_H