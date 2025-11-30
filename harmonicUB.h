#ifndef HARMONICUB_H
#define HARMONICUB_H


#include "atom.h"


class HarmonicUB
{
    public:

    HarmonicUB(float UB_force_constant, float UB_force_equil) : 
    UB_force_constant(UB_force_constant), 
    UB_force_equil(UB_force_equil) {}
    
    HarmonicUB(float UB_force_constant, float UB_force_equil, Atom atomA, Atom atomB) : 
    UB_force_constant(UB_force_constant), 
    UB_force_equil(UB_force_equil), 
    atomA(atomA), 
    atomB(atomB) {}

    HarmonicUB(Atom atomA, Atom atomB) : 
    UB_force_constant(0.0), 
    UB_force_equil(0.0), 
    atomA(atomA), 
    atomB(atomB) {}

    float return_force(float distance);


    private:
        float UB_force_constant;
        float UB_force_equil;
        Atom atomA;
        Atom atomB;
};
#endif //HARMONICUB_H