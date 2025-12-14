#ifndef HARMONIIMPROPER_H
#define HARMONIIMPROPER_H


// #include "atom.h"


class HarmonicImproper
{

    public:
        HarmonicImproper(const double IMP_force_constant, const double IMP_Phase) :
            IMP_force_constant(IMP_force_constant),
            IMP_Phase(IMP_Phase),
            atomA_index(-1),
            atomB_index(-1),
            atomC_index(-1),
            atomD_index(-1), type(0), dihedral(0), energy(0) {}

        HarmonicImproper(const double IMP_force_constant, const double IMP_Phase, const int atomA_index, const int atomB_index, const int atomC_index, const int atomD_index) :
            IMP_force_constant(IMP_force_constant),
            IMP_Phase(IMP_Phase),
            atomA_index(atomA_index),
            atomB_index(atomB_index),
            atomC_index(atomC_index),
            atomD_index(atomD_index), type(0), dihedral(0), energy(0) {}

        HarmonicImproper(const int atomA_index, const int atomB_index, const int atomC_index, const int atomD_index) :
            IMP_force_constant(0.0),
            IMP_Phase(0.0),
            atomA_index(atomA_index),
            atomB_index(atomB_index),
            atomC_index(atomC_index),
            atomD_index(atomD_index), type(0),
            dihedral(0), energy(0) {}

        [[nodiscard]] double return_energy(double distance) const;
        void set_energy(const double energy_) {this->energy = energy_;}


        void set_type(const int type_id) { this->type = type_id; }
        [[nodiscard]] int get_type() const { return type; }

        void set_IMP_force_constant(const double force_constant) {this->IMP_force_constant = force_constant; }
        [[nodiscard]] double get_IMP_force_constant() const {return IMP_force_constant; }

        void set_IMP_phase_value(const double phase_value) {this->IMP_Phase = phase_value; }
        [[nodiscard]] double get_IMP_phase_value() const {return IMP_Phase; }

        void set_imp_dihedral(const double dihedral_) {this->dihedral = dihedral_;}
        [[nodiscard]] double get_imp_dihedral() const {return this->dihedral;}

        [[nodiscard]] int get_atomA_index() const {return atomA_index; }
        [[nodiscard]] int get_atomB_index() const {return atomB_index; }
        [[nodiscard]] int get_atomC_index() const {return atomC_index; }
        [[nodiscard]] int get_atomD_index() const {return atomD_index; }
    
    private:
        double IMP_force_constant;
        double IMP_Phase;
        int atomA_index;
        int atomB_index;
        int atomC_index;
        int atomD_index;
        int type;
        double dihedral;
        double energy;
};
#endif //HARMONIIMPROPER_H