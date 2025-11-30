#ifndef HARMONIIMPROPER_H
#define HARMONIIMPROPER_H


#include "atom.h"


class HarmonicImproper
{

    public:
        HarmonicImproper(float IMP_force_constant, float IMP_force_equil) : 
        IMP_force_constant(IMP_force_constant), 
        IMP_force_equil(IMP_force_equil), type(0) {}
        
        HarmonicImproper(float IMP_force_constant, float IMP_force_equil, Atom atomA, Atom atomB, Atom atomC, Atom atomD) : 
        IMP_force_constant(IMP_force_constant), 
        IMP_force_equil(IMP_force_equil), 
        atomA(atomA), 
        atomB(atomB),
        atomC(atomC),
        atomD(atomD), type(0) {}

        HarmonicImproper(Atom atomA, Atom atomB, Atom atomC, Atom atomD) : 
        IMP_force_constant(0.0), 
        IMP_force_equil(0.0), 
        atomA(atomA), 
        atomB(atomB), 
        atomC(atomC),
        atomD(atomD), type(0) {}

        float return_force(float distance);

        void set_type(int type_id) { type = type_id; }
        int get_type() const { return type; }

        void set_IMP_force_constant(float force_constant) {IMP_force_constant = force_constant; }
        float get_IMP_force_constant() const {return IMP_force_constant; }

        void set_IMP_equil_value(float equil_value) {IMP_force_equil = equil_value; }
        float get_IMP_equil_value() const {return IMP_force_equil; }

        const Atom& get_atomA() const {return atomA; }
        const Atom& get_atomB() const {return atomB; }
        const Atom& get_atomC() const {return atomC; }
        const Atom& get_atomD() const {return atomD; }
    
    private:
        float IMP_force_constant;
        float IMP_force_equil;
        Atom atomA;
        Atom atomB;
        Atom atomC;
        Atom atomD;
        int type;
};
#endif //HARMONIIMPROPER_H