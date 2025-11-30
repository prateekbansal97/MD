#ifndef HARMONICUB_H
#define HARMONICUB_H


#include "atom.h"


class HarmonicUB
{

    public:
        HarmonicUB(float UB_force_constant, float UB_force_equil) : 
        UB_force_constant(UB_force_constant), 
        UB_force_equil(UB_force_equil), type(0) {}
        
        HarmonicUB(float UB_force_constant, float UB_force_equil, Atom atomA, Atom atomB) : 
        UB_force_constant(UB_force_constant), 
        UB_force_equil(UB_force_equil), 
        atomA(atomA), 
        atomB(atomB), type(0) {}

        HarmonicUB(Atom atomA, Atom atomB) : 
        UB_force_constant(0.0), 
        UB_force_equil(0.0), 
        atomA(atomA), 
        atomB(atomB), type(0) {}

        float return_force(float distance);

        void set_type(int type_id) { type = type_id; }
        int get_type() const { return type; }

        void set_UB_force_constant(float force_constant) {UB_force_constant = force_constant; }
        float get_UB_force_constant() const {return UB_force_constant; }

        void set_UB_equil_value(float equil_value) {UB_force_equil = equil_value; }
        float get_UB_equil_value() const {return UB_force_equil; }

        const Atom& get_atomA() const {return atomA; }
        const Atom& get_atomB() const {return atomB; }
    
        private:
        float UB_force_constant;
        float UB_force_equil;
        Atom atomA;
        Atom atomB;
        int type;
};
#endif //HARMONICUB_H