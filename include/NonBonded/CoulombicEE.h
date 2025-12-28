//
// Created by Prateek Bansal on 12/13/25.
//

#ifndef MD_COULOMBICEE_H
#define MD_COULOMBICEE_H

namespace md { class System; }

namespace md::NonBonded
{
    class CoulombicEE
    {
    private:
        static double CalculateEnergy(double r, double chargeA, double chargeB, double epsilon);
        static double CalculateGradient(double r, double chargeA, double chargeB, double epsilon);
        static double CalculateEnergyCutoff(double r, double rcutoff, double chargeA, double chargeB, double epsilon);
        static double CalculateEnergyEwaldDirectTerm(double r, double chargeA, double chargeB, double ewald_alpha);
        static double CalculateEnergyEwaldDirectTerm_14(double r, double chargeA, double chargeB, double ewald_alpha, double scale);



        friend class md::System;
    };
}
#endif //MD_COULOMBICEE_H