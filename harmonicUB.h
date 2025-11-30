#ifndef HARMONICUB_H
#define HARMONICUB_H


#include "atom.h"


class HarmonicUB
{

    public:
        HarmonicUB(double UB_force_constant, double UB_force_equil) : 
        UB_force_constant(UB_force_constant), 
        UB_force_equil(UB_force_equil), type(0) {}
        
        HarmonicUB(double UB_force_constant, double UB_force_equil, Atom atomA, Atom atomB) : 
        UB_force_constant(UB_force_constant), 
        UB_force_equil(UB_force_equil), 
        atomA(atomA), 
        atomB(atomB), type(0) {}

        HarmonicUB(Atom atomA, Atom atomB) : 
        UB_force_constant(0.0), 
        UB_force_equil(0.0), 
        atomA(atomA), 
        atomB(atomB), type(0) {}

        double return_energy(double distance);

        void set_type(int type_id) { type = type_id; }
        int get_type() const { return type; }

        void set_UB_force_constant(double force_constant) {UB_force_constant = force_constant; }
        double get_UB_force_constant() const {return UB_force_constant; }

        void set_UB_equil_value(double equil_value) {UB_force_equil = equil_value; }
        double get_UB_equil_value() const {return UB_force_equil; }

        const Atom& get_atomA() const {return atomA; }
        const Atom& get_atomB() const {return atomB; }
    
        private:
        double UB_force_constant;
        double UB_force_equil;
        Atom atomA;
        Atom atomB;
        int type;
};
#endif //HARMONICUB_H