#ifndef HARMONICUB_H
#define HARMONICUB_H




class HarmonicUB
{

    public:
        HarmonicUB(const double UB_force_constant, const double UB_force_equil) :
            UB_force_constant_(UB_force_constant),
            UB_force_equil_(UB_force_equil),
            atomA_index_(-1),
            atomB_index_(-1),
            type_(0), distance_(0), energy_(0) {}

        HarmonicUB(const double UB_force_constant, const double UB_force_equil, const int atomA_index, const int atomB_index) :
            UB_force_constant_(UB_force_constant),
            UB_force_equil_(UB_force_equil),
            atomA_index_(atomA_index),
            atomB_index_(atomB_index),
            type_(0), distance_(0), energy_(0) {}

        HarmonicUB(const int atomA_index, const int atomB_index) :
            UB_force_constant_(0.0),
            UB_force_equil_(0.0),
            atomA_index_(atomA_index),
            atomB_index_(atomB_index),
            type_(0), distance_(0), energy_(0) {}

        // double return_energy(double distance);
        void set_energy(const double energy) {this->energy_ = energy;}
        void set_type(const int type_id) { this->type_ = type_id; }
        void set_UB_force_constant(const double force_constant) {this->UB_force_constant_ = force_constant; }
        void set_UB_equil_value(const double equil_value) {this->UB_force_equil_ = equil_value; }
        void set_distance_value(const double distance) {this->distance_ = distance;}

        [[nodiscard]] int get_type() const { return type_; }
        [[nodiscard]] double get_UB_force_constant() const {return UB_force_constant_; }
        [[nodiscard]] double get_UB_equil_value() const {return UB_force_equil_; }
        [[nodiscard]] int get_atomA_index() const {return atomA_index_; }
        [[nodiscard]] int get_atomB_index() const {return atomB_index_; }

    private:
        double UB_force_constant_;
        double UB_force_equil_;
        int atomA_index_;
        int atomB_index_;
        int type_;
        double distance_;
        double energy_;

        [[nodiscard]] double calculate_energy(double distance) const;

        friend class System;

};
#endif //HARMONICUB_H