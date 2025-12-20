#ifndef HARMONICANGLE_H
#define HARMONICANGLE_H

class HarmonicAngle
{

    public:
        HarmonicAngle(const double Angle_force_constant, const double Angle_Equil) :
            Angle_force_constant_(Angle_force_constant),
            Angle_Equil_(Angle_Equil), atomA_index_(-1),
            atomB_index_(-1),
            atomC_index_(-1),
            type_(0),
            isH_(false), angle_(0), energy_(0) {}

        HarmonicAngle(const double Angle_force_constant, const double Angle_Equil, const int atomA_index, const int atomB_index, const int atomC_index, const bool isH) :
            Angle_force_constant_(Angle_force_constant),
            Angle_Equil_(Angle_Equil),
            atomA_index_(atomA_index),
            atomB_index_(atomB_index), atomC_index_(atomC_index),
            type_(0), isH_(isH), angle_(0), energy_(0) {}

        HarmonicAngle(const int atomA_index, const int atomB_index, const int atomC_index) :
            Angle_force_constant_(0.0),
            Angle_Equil_(0.0),
            atomA_index_(atomA_index),
            atomB_index_(atomB_index), atomC_index_(atomC_index), type_(0), isH_(false),
            angle_(0), energy_(0) {}

        void set_energy(const double energy_) {this->energy_ = energy_;}
        void set_type(const int type_id) { this->type_ = type_id; }
        void set_Angle_force_constant(const double force_constant) {this->Angle_force_constant_ = force_constant; }
        void set_Angle_equil_angle(const double equil_angle) {this->Angle_Equil_ = equil_angle; }
        void set_angle(const double angle_) { this->angle_ = angle_;}

        [[nodiscard]] double get_Angle_equil_angle() const {return Angle_Equil_; }
        [[nodiscard]] double get_Angle_force_constant() const {return Angle_force_constant_; }
        [[nodiscard]] int get_type() const { return type_; }
        [[nodiscard]] int get_atomA_index() const {return atomA_index_; }
        [[nodiscard]] int get_atomB_index() const {return atomB_index_; }
        [[nodiscard]] int get_atomC_index() const {return atomC_index_; }
        [[nodiscard]] bool get_isH() const {return isH_;}
        [[nodiscard]] double get_angle() const {return angle_;}
        [[nodiscard]] double get_energy() const {return energy_;}


    private:
        [[nodiscard]] double calculate_energy(double angle) const;

        double Angle_force_constant_;
        double Angle_Equil_;
        int atomA_index_;
        int atomB_index_;
        int atomC_index_;
        int type_;
        bool isH_;
        double angle_;
        double energy_;
        friend class System;
};
#endif //HARMONICANGLE_H