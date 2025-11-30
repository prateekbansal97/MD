#ifndef HARMONIIMPROPER_H
#define HARMONIIMPROPER_H


#include "atom.h"


class HarmonicImproper
{

    public:
        HarmonicImproper(double IMP_force_constant, double IMP_Phase) : 
        IMP_force_constant(IMP_force_constant), 
        IMP_Phase(IMP_Phase), type(0) {}
        
        HarmonicImproper(double IMP_force_constant, double IMP_Phase, Atom atomA, Atom atomB, Atom atomC, Atom atomD) : 
        IMP_force_constant(IMP_force_constant), 
        IMP_Phase(IMP_Phase), 
        atomA(atomA), 
        atomB(atomB),
        atomC(atomC),
        atomD(atomD), type(0) {}

        HarmonicImproper(Atom atomA, Atom atomB, Atom atomC, Atom atomD) : 
        IMP_force_constant(0.0), 
        IMP_Phase(0.0), 
        atomA(atomA), 
        atomB(atomB), 
        atomC(atomC),
        atomD(atomD), type(0) {}

        double return_energy(double distance);

        void set_type(int type_id) { type = type_id; }
        int get_type() const { return type; }

        void set_IMP_force_constant(double force_constant) {IMP_force_constant = force_constant; }
        double get_IMP_force_constant() const {return IMP_force_constant; }

        void set_IMP_phase_value(double phase_value) {IMP_Phase = phase_value; }
        double get_IMP_phase_value() const {return IMP_Phase; }

        const Atom& get_atomA() const {return atomA; }
        const Atom& get_atomB() const {return atomB; }
        const Atom& get_atomC() const {return atomC; }
        const Atom& get_atomD() const {return atomD; }
    
    private:
        double IMP_force_constant;
        double IMP_Phase;
        Atom atomA;
        Atom atomB;
        Atom atomC;
        Atom atomD;
        int type;
};
#endif //HARMONIIMPROPER_H