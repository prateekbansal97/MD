#ifndef COSINEDIHEDRAL_H
#define COSINEDIHEDRAL_H

namespace md { class System; }

namespace md::Bonded
{
    class CosineDihedral
    {

    public:
        CosineDihedral(const double Dihedral_force_constant, const double Dihedral_Phase, const double Dihedral_Periodicity,
                        const bool improper = false, const bool exclude_14 = false) :
        Dihedral_force_constant_(Dihedral_force_constant),
        Dihedral_Phase_(Dihedral_Phase),
        Dihedral_Periodicity_(Dihedral_Periodicity),
        atomA_index_(-1),
        atomB_index_(-1),
        atomC_index_(-1),
        atomD_index_(-1),
        type_(0),
        isH_(false),
        improper_(improper),
        exclude_14_(exclude_14) {}

        CosineDihedral(const double Dihedral_force_constant, const double Dihedral_Phase, const double Dihedral_Periodicity,
                        const int atomA_index, const int atomB_index, const int atomC_index, const int atomD_index,
                        const bool isH, const bool improper = false, const bool exclude_14 = false) :
        Dihedral_force_constant_(Dihedral_force_constant),
        Dihedral_Phase_(Dihedral_Phase),
        Dihedral_Periodicity_(Dihedral_Periodicity),
        atomA_index_(atomA_index),
        atomB_index_(atomB_index),
        atomC_index_(atomC_index),
        atomD_index_(atomD_index),
        type_(0),
        isH_(isH),
        improper_(improper),
        exclude_14_(exclude_14) {}

        CosineDihedral(const int atomA_index, const int atomB_index, const int atomC_index, const int atomD_index, const bool improper = false, const bool exclude_14 = false) :
        Dihedral_force_constant_(0.0),
        Dihedral_Phase_(0.0),
        Dihedral_Periodicity_(0.0),
        atomA_index_(atomA_index),
        atomB_index_(atomB_index),
        atomC_index_(atomC_index),
        atomD_index_(atomD_index),
        type_(0),
        isH_(false),
        improper_(improper),
        exclude_14_(exclude_14) {}

        void set_energy(const double energ) {this->energy_ = energ;}
        void set_type(const int type_id) { this->type_ = type_id; }
        void set_Dihedral_force_constant(const double force_constant) {this->Dihedral_force_constant_ = force_constant; }
        void set_Dihedral_phase(const double phase) {this->Dihedral_Phase_ = phase; }
        void set_Dihedral_Periodicity(const double Periodicity) {this->Dihedral_Periodicity_ = Periodicity; }
        void set_cosine_dihedral(const double dihedra) {this->dihedral_ = dihedra;}

        [[nodiscard]] int get_type() const { return type_; }
        [[nodiscard]] double get_Dihedral_force_constant() const {return Dihedral_force_constant_; }
        [[nodiscard]] double get_Dihedral_phase() const {return Dihedral_Phase_; }
        [[nodiscard]] double get_Dihedral_Periodicity() const {return Dihedral_Periodicity_; }
        [[nodiscard]] int get_atomA_index() const {return atomA_index_; }
        [[nodiscard]] int get_atomB_index() const {return atomB_index_; }
        [[nodiscard]] int get_atomC_index() const {return atomC_index_; }
        [[nodiscard]] int get_atomD_index() const {return atomD_index_; }
        [[nodiscard]] bool get_isH() const {return isH_;}
        [[nodiscard]] bool get_improper() const {return improper_;}
        [[nodiscard]] bool get_exclude_14() const {return exclude_14_;}
        [[nodiscard]] double get_energy() const {return energy_; }
        [[nodiscard]] double get_dihedral() const {return dihedral_; }


    private:
        [[nodiscard]] double calculate_energy(double dihedra) const;
        double Dihedral_force_constant_;
        double Dihedral_Phase_;
        double Dihedral_Periodicity_;
        int atomA_index_;
        int atomB_index_;
        int atomC_index_;
        int atomD_index_;
        int type_;
        bool isH_;
        bool improper_;
        bool exclude_14_;
        double dihedral_{};
        double energy_{};

        friend class ::md::System;
    };
}
#endif //COSINEDIHEDRAL_H