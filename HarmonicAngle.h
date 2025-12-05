#ifndef HARMONICANGLE_H
#define HARMONICANGLE_H

class HarmonicAngle
{

    public:
        HarmonicAngle(double Angle_force_constant, double Angle_Equil) : 
        Angle_force_constant(Angle_force_constant), 
        Angle_Equil(Angle_Equil), type(0), isH(false) {}
        
        HarmonicAngle(double Angle_force_constant, double Angle_Equil, int atomA_index, int atomB_index, int atomC_index, bool isH) : 
        Angle_force_constant(Angle_force_constant), 
        Angle_Equil(Angle_Equil), 
        atomA_index(atomA_index), 
        atomB_index(atomB_index), atomC_index(atomC_index), type(0), isH(isH) {}

        HarmonicAngle(int atomA_index, int atomB_index, int atomC_index) : 
        Angle_force_constant(0.0), 
        Angle_Equil(0.0), 
        atomA_index(atomA_index), 
        atomB_index(atomB_index), atomC_index(atomC_index), type(0), isH(false) {}

        double return_energy(double angle);

        void set_type(int type_id) { type = type_id; }
        int get_type() const { return type; }

        void set_Angle_force_constant(double force_constant) {Angle_force_constant = force_constant; }
        double get_Angle_force_constant() const {return Angle_force_constant; }

        void set_Angle_equil_angle(double equil_angle) {Angle_Equil = equil_angle; }
        double get_Angle_equil_angle() const {return Angle_Equil; }

        const int get_atomA_index() const {return atomA_index; }
        const int get_atomB_index() const {return atomB_index; }
        const int get_atomC_index() const {return atomC_index; }

        const bool get_isH() const {return isH;}
    
    private:
        double Angle_force_constant;
        double Angle_Equil;
        int atomA_index;
        int atomB_index;
        int atomC_index;
        int type;
        bool isH;
};
#endif //HARMONICANGLE_H