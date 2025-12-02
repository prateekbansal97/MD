#ifndef COSINEDIHEDRAL_H
#define COSINEDIHEDRAL_H


#include "atom.h"


class CosineDihedral
{

    public:
        CosineDihedral(double Dihedral_force_constant, double Dihedral_Phase, double Dihedral_Periodicity, 
                        bool improper = false, bool exclude_14 = false) : 
        Dihedral_force_constant(Dihedral_force_constant), 
        Dihedral_Phase(Dihedral_Phase),
        Dihedral_Periodicity(Dihedral_Periodicity), 
        type(0), 
        improper(improper),
        exclude_14(exclude_14),
        isH(false) {}
        
        CosineDihedral(double Dihedral_force_constant, double Dihedral_Phase, double Dihedral_Periodicity, 
                        Atom atomA, Atom atomB, Atom atomC, Atom atomD, 
                        bool isH, bool improper = false, bool exclude_14 = false) : 
        Dihedral_force_constant(Dihedral_force_constant), 
        Dihedral_Phase(Dihedral_Phase),
        Dihedral_Periodicity(Dihedral_Periodicity),
        atomA(atomA), 
        atomB(atomB), 
        atomC(atomC),
        atomD(atomD), 
        type(0), 
        improper(improper),
        exclude_14(exclude_14),
        isH(isH) {}

        CosineDihedral(Atom atomA, Atom atomB, Atom atomC, Atom atomD, bool improper = false, bool exclude_14 = false) : 
        Dihedral_force_constant(0.0), 
        Dihedral_Phase(0.0),
        Dihedral_Periodicity(0.0), 
        atomA(atomA), 
        atomB(atomB), 
        atomC(atomC), 
        atomD(atomD),
        type(0), 
        improper(improper),
        exclude_14(exclude_14),
        isH(false) {}

        double return_energy(double dihedral);

        void set_type(int type_id) { type = type_id; }
        int get_type() const { return type; }

        void set_Dihedral_force_constant(double force_constant) {Dihedral_force_constant = force_constant; }
        double get_Dihedral_force_constant() const {return Dihedral_force_constant; }

        void set_Dihedral_phase(double phase) {Dihedral_Phase = phase; }
        double get_Dihedral_phase() const {return Dihedral_Phase; }

        void set_Dihedral_Periodicity(double Periodicity) {Dihedral_Periodicity = Periodicity; }
        double get_Dihedral_Periodicity() const {return Dihedral_Periodicity; }

        const Atom& get_atomA() const {return atomA; }
        const Atom& get_atomB() const {return atomB; }
        const Atom& get_atomC() const {return atomC; }
        const Atom& get_atomD() const {return atomD; }

        const bool get_isH() const {return isH;}
    
    private:
        double Dihedral_force_constant;
        double Dihedral_Phase;
        double Dihedral_Periodicity;
        Atom atomA;
        Atom atomB;
        Atom atomC;
        Atom atomD;
        int type;
        bool isH;
        bool improper;
        bool exclude_14;
};
#endif //COSINEDIHEDRAL_H