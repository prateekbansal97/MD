#ifndef HARMONICBOND_H
#define HARMONICBOND_H

namespace md { class System; }

namespace md::Bonded
{
    class HarmonicBond
    {

    public:
        HarmonicBond(const double Bond_force_constant, const double Bond_Equil) :
            Bond_force_constant_(Bond_force_constant),
            Bond_Equil_(Bond_Equil), atomA_index_(-1),
            atomB_index_(-1),
            type_(0),
            isH_(false), distance_(0), energy_(0) {}

        HarmonicBond(const double Bond_force_constant, const double Bond_Equil, const int atomA_index, const int atomB_index, const bool isH) :
            Bond_force_constant_(Bond_force_constant),
            Bond_Equil_(Bond_Equil),
            atomA_index_(atomA_index),
            atomB_index_(atomB_index), type_(0),
            isH_(isH), distance_(0), energy_(0) {}

        HarmonicBond(const int atomA_index, const int atomB_index) :
            Bond_force_constant_(0.0),
            Bond_Equil_(0.0),
            atomA_index_(atomA_index),
            atomB_index_(atomB_index), type_(0),
            isH_(false), distance_(0), energy_(0) {}

        void set_energy(const double energy_) { this->energy_ = energy_;}
        void set_type(const int type_id) { this->type_ = type_id; }
        void set_Bond_force_constant(const double force_constant) {this->Bond_force_constant_ = force_constant; }
        void set_Bond_equil_length(const double equil_length) {this->Bond_Equil_ = equil_length; }
        void set_distance(const double distance_) {this->distance_ = distance_;}

        [[nodiscard]] int get_type() const { return type_; }
        [[nodiscard]] double get_Bond_force_constant() const {return Bond_force_constant_; }
        [[nodiscard]] double get_Bond_equil_length() const {return Bond_Equil_; }
        [[nodiscard]] int get_atomA_index() const {return atomA_index_; }
        [[nodiscard]] int get_atomB_index() const {return atomB_index_; }
        [[nodiscard]] bool get_isH() const {return isH_;}
        [[nodiscard]] double get_distance() const {return distance_;}
        [[nodiscard]] double get_energy() const {return energy_;}



    private:
        [[nodiscard]] double calculate_energy(double distance) const;

        double Bond_force_constant_;
        double Bond_Equil_;
        int atomA_index_;
        int atomB_index_;
        int type_;
        bool isH_;
        double distance_;
        double energy_;

        friend class ::md::System;
    };
} //namespace md
#endif //HARMONICBOND_H