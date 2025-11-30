#ifndef HARMONIBOND_H
#define HARMONIBOND_H


#include "atom.h"


class HarmonicBond
{

    public:
        HarmonicBond(double Bond_force_constant, double Bond_Equil) : 
        Bond_force_constant(Bond_force_constant), 
        Bond_Equil(Bond_Equil), type(0) {}
        
        HarmonicBond(double Bond_force_constant, double Bond_Equil, Atom atomA, Atom atomB) : 
        Bond_force_constant(Bond_force_constant), 
        Bond_Equil(Bond_Equil), 
        atomA(atomA), 
        atomB(atomB), type(0) {}

        HarmonicBond(Atom atomA, Atom atomB) : 
        Bond_force_constant(0.0), 
        Bond_Equil(0.0), 
        atomA(atomA), 
        atomB(atomB), type(0) {}

        double return_energy(double distance);

        void set_type(int type_id) { type = type_id; }
        int get_type() const { return type; }

        void set_Bond_force_constant(double force_constant) {Bond_force_constant = force_constant; }
        double get_Bond_force_constant() const {return Bond_force_constant; }

        void set_Bond_phase_value(double phase_value) {Bond_Equil = phase_value; }
        double get_Bond_phase_value() const {return Bond_Equil; }

        const Atom& get_atomA() const {return atomA; }
        const Atom& get_atomB() const {return atomB; }
    
    private:
        double Bond_force_constant;
        double Bond_Equil;
        Atom atomA;
        Atom atomB;
        int type;
};
#endif //HARMONIBOND_H