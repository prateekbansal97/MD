//
// Created by Prateek Bansal on 12/12/25.
//

#ifndef MD_LENNARDJONES_H
#define MD_LENNARDJONES_H

class LennardJones
{
    public:
    LennardJones();
    ~LennardJones();

    static double CalculateEnergy(double distance, double Aij, double Bij);
    static double CalculateGradient(double r2, double Aij, double Bij);

};

#endif //MD_LENNARDJONES_H