#ifndef HARMONICUB_H
#define HARMONICUB_H


// #include "atom.h"


class HarmonicUB
{

    public:
        HarmonicUB(double UB_force_constant, double UB_force_equil) : 
        UB_force_constant(UB_force_constant), 
        UB_force_equil(UB_force_equil),
        atomA_index(-1), 
        atomB_index(-1), type(0) {}
        
        HarmonicUB(double UB_force_constant, double UB_force_equil, int atomA_index, int atomB_index) : 
        UB_force_constant(UB_force_constant), 
        UB_force_equil(UB_force_equil), 
        atomA_index(atomA_index), 
        atomB_index(atomB_index), type(0) {}

        HarmonicUB(int atomA_index, int atomB_index) : 
        UB_force_constant(0.0), 
        UB_force_equil(0.0), 
        atomA_index(atomA_index), 
        atomB_index(atomB_index), type(0) {}

        double return_energy(double distance);
        void set_energy(double energy) {this->energy = energy;}

        void set_type(int type_id) { this->type = type_id; }
        int get_type() const { return type; }

        void set_UB_force_constant(double force_constant) {this->UB_force_constant = force_constant; }
        double get_UB_force_constant() const {return UB_force_constant; }

        void set_UB_equil_value(double equil_value) {this->UB_force_equil = equil_value; }
        double get_UB_equil_value() const {return UB_force_equil; }

        void set_distance_value(double distance) {this->distance = distance;}

        const int get_atomA_index() const {return atomA_index; }
        const int get_atomB_index() const {return atomB_index; }
    
        private:
        double UB_force_constant;
        double UB_force_equil;
        int atomA_index;
        int atomB_index;
        int type;
        double distance;
        double energy;
};
#endif //HARMONICUB_H