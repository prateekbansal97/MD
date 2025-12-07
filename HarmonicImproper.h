#ifndef HARMONIIMPROPER_H
#define HARMONIIMPROPER_H


// #include "atom.h"


class HarmonicImproper
{

    public:
        HarmonicImproper(double IMP_force_constant, double IMP_Phase) : 
        IMP_force_constant(IMP_force_constant), 
        IMP_Phase(IMP_Phase), 
        atomA_index(-1), 
        atomB_index(-1), 
        atomC_index(-1),
        atomD_index(-1), type(0) {}
        
        HarmonicImproper(double IMP_force_constant, double IMP_Phase, int atomA_index, int atomB_index, int atomC_index, int atomD_index) : 
        IMP_force_constant(IMP_force_constant), 
        IMP_Phase(IMP_Phase), 
        atomA_index(atomA_index), 
        atomB_index(atomB_index), 
        atomC_index(atomC_index),
        atomD_index(atomD_index), type(0) {}

        HarmonicImproper(int atomA_index, int atomB_index, int atomC_index, int atomD_index) : 
        IMP_force_constant(0.0), 
        IMP_Phase(0.0), 
        atomA_index(atomA_index), 
        atomB_index(atomB_index), 
        atomC_index(atomC_index),
        atomD_index(atomD_index), type(0) {}

        double return_energy(double distance);

        void set_type(int type_id) { this->type = type_id; }
        int get_type() const { return type; }

        void set_IMP_force_constant(double force_constant) {this->IMP_force_constant = force_constant; }
        double get_IMP_force_constant() const {return IMP_force_constant; }

        void set_IMP_phase_value(double phase_value) {this->IMP_Phase = phase_value; }
        double get_IMP_phase_value() const {return IMP_Phase; }

        const int get_atomA_index() const {return atomA_index; }
        const int get_atomB_index() const {return atomB_index; }
        const int get_atomC_index() const {return atomC_index; }
        const int get_atomD_index() const {return atomD_index; }
    
    private:
        double IMP_force_constant;
        double IMP_Phase;
        int atomA_index;
        int atomB_index;
        int atomC_index;
        int atomD_index;
        int type;
};
#endif //HARMONIIMPROPER_H