#ifndef HARMONIANGLE_H
#define HARMONIANGLE_H


#include "atom.h"


class HarmonicAngle
{

    public:
        HarmonicAngle(double Angle_force_constant, double Angle_Equil) : 
        Angle_force_constant(Angle_force_constant), 
        Angle_Equil(Angle_Equil), type(0), isH(false) {}
        
        HarmonicAngle(double Angle_force_constant, double Angle_Equil, Atom atomA, Atom atomB, Atom atomC, bool isH) : 
        Angle_force_constant(Angle_force_constant), 
        Angle_Equil(Angle_Equil), 
        atomA(atomA), 
        atomB(atomB), atomC(atomC), type(0), isH(isH) {}

        HarmonicAngle(Atom atomA, Atom atomB, Atom atomC) : 
        Angle_force_constant(0.0), 
        Angle_Equil(0.0), 
        atomA(atomA), 
        atomB(atomB), atomC(atomC), type(0), isH(false) {}

        double return_energy(double angle);

        void set_type(int type_id) { type = type_id; }
        int get_type() const { return type; }

        void set_Angle_force_constant(double force_constant) {Angle_force_constant = force_constant; }
        double get_Angle_force_constant() const {return Angle_force_constant; }

        void set_Angle_equil_angle(double equil_angle) {Angle_Equil = equil_angle; }
        double get_Angle_equil_angle() const {return Angle_Equil; }

        const Atom& get_atomA() const {return atomA; }
        const Atom& get_atomB() const {return atomB; }
        const Atom& get_atomC() const {return atomC; }

        const bool get_isH() const {return isH;}
    
    private:
        double Angle_force_constant;
        double Angle_Equil;
        Atom atomA;
        Atom atomB;
        Atom atomC;
        int type;
        bool isH;
};
#endif //HARMONIANGLE_H