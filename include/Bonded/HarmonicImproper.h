#ifndef HARMONIIMPROPER_H
#define HARMONIIMPROPER_H

namespace md { class System; }

namespace md::Bonded
{
    class HarmonicImproper
    {

    public:
        HarmonicImproper(const double IMP_force_constant, const double IMP_Phase) :
            IMP_force_constant_(IMP_force_constant),
            IMP_Phase_(IMP_Phase),
            atomA_index_(-1),
            atomB_index_(-1),
            atomC_index_(-1),
            atomD_index_(-1), type_(0), dihedral_(0), energy_(0) {}

        HarmonicImproper(const double IMP_force_constant, const double IMP_Phase, const int atomA_index, const int atomB_index, const int atomC_index, const int atomD_index) :
            IMP_force_constant_(IMP_force_constant),
            IMP_Phase_(IMP_Phase),
            atomA_index_(atomA_index),
            atomB_index_(atomB_index),
            atomC_index_(atomC_index),
            atomD_index_(atomD_index), type_(0), dihedral_(0), energy_(0) {}

        HarmonicImproper(const int atomA_index, const int atomB_index, const int atomC_index, const int atomD_index) :
            IMP_force_constant_(0.0),
            IMP_Phase_(0.0),
            atomA_index_(atomA_index),
            atomB_index_(atomB_index),
            atomC_index_(atomC_index),
            atomD_index_(atomD_index), type_(0),
            dihedral_(0), energy_(0) {}

        void set_energy(const double energy) {this->energy_ = energy;}
        void set_type(const int type_id) { this->type_ = type_id; }
        void set_IMP_force_constant(const double force_constant) {this->IMP_force_constant_ = force_constant; }
        void set_IMP_phase_value(const double phase_value) {this->IMP_Phase_ = phase_value; }
        void set_imp_dihedral(const double dihedral) {this->dihedral_ = dihedral;}

        [[nodiscard]] int get_type() const { return type_; }
        [[nodiscard]] double get_IMP_force_constant() const {return IMP_force_constant_; }
        [[nodiscard]] double get_IMP_phase_value() const {return IMP_Phase_; }
        [[nodiscard]] double get_imp_dihedral() const {return this->dihedral_;}
        [[nodiscard]] int get_atomA_index() const {return atomA_index_; }
        [[nodiscard]] int get_atomB_index() const {return atomB_index_; }
        [[nodiscard]] int get_atomC_index() const {return atomC_index_; }
        [[nodiscard]] int get_atomD_index() const {return atomD_index_; }
        [[nodiscard]] double get_energy() const {return energy_;}

    private:
        double IMP_force_constant_;
        double IMP_Phase_;
        int atomA_index_;
        int atomB_index_;
        int atomC_index_;
        int atomD_index_;
        int type_;
        double dihedral_;
        double energy_;

        [[nodiscard]] double calculate_energy(double distance) const;

        friend class ::md::System;

    };
} //namespace md
#endif //HARMONIIMPROPER_H