#ifndef HARMONICBOND_H
#define HARMONICBOND_H


// #include "atom.h"


class HarmonicBond
{

    public:
        HarmonicBond(double Bond_force_constant, double Bond_Equil) : 
        Bond_force_constant(Bond_force_constant), 
        Bond_Equil(Bond_Equil), type(0), 
        atomA_index(-1), 
        atomB_index(-1), 
        isH(false) {}
        
        HarmonicBond(double Bond_force_constant, double Bond_Equil, int atomA_index, int atomB_index, bool isH) : 
        Bond_force_constant(Bond_force_constant), 
        Bond_Equil(Bond_Equil), 
        atomA_index(atomA_index), 
        atomB_index(atomB_index), type(0), isH(isH) {}

        HarmonicBond(int atomA_index, int atomB_index) : 
        Bond_force_constant(0.0), 
        Bond_Equil(0.0), 
        atomA_index(atomA_index), 
        atomB_index(atomB_index), type(0), isH(false) {}

        double return_energy(double distance);

        void set_type(int type_id) { this->type = type_id; }
        int get_type() const { return type; }

        void set_Bond_force_constant(double force_constant) {this->Bond_force_constant = force_constant; }
        double get_Bond_force_constant() const {return Bond_force_constant; }

        void set_Bond_equil_length(double equil_length) {this->Bond_Equil = equil_length; }
        double get_Bond_equil_length() const {return Bond_Equil; }

        void set_distance(double distance) {this->distance = distance;}

        const int get_atomA_index() const {return atomA_index; }
        const int get_atomB_index() const {return atomB_index; }

        const bool get_isH() const {return isH;}


    private:
        double Bond_force_constant;
        double Bond_Equil;
        int atomA_index;
        int atomB_index;
        int type;
        bool isH;
        double distance;
};
#endif //HARMONICBOND_H