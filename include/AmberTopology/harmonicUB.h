#ifndef HARMONICUB_H
#define HARMONICUB_H




class HarmonicUB
{

    public:
        HarmonicUB(const double UB_force_constant, const double UB_force_equil) :
            UB_force_constant(UB_force_constant),
            UB_force_equil(UB_force_equil),
            atomA_index(-1),
            atomB_index(-1),
            type(0), distance(0), energy(0) {}

        HarmonicUB(const double UB_force_constant, const double UB_force_equil, const int atomA_index, const int atomB_index) :
            UB_force_constant(UB_force_constant),
            UB_force_equil(UB_force_equil),
            atomA_index(atomA_index),
            atomB_index(atomB_index),
            type(0), distance(0), energy(0) {}

        HarmonicUB(const int atomA_index, const int atomB_index) :
            UB_force_constant(0.0),
            UB_force_equil(0.0),
            atomA_index(atomA_index),
            atomB_index(atomB_index),
            type(0), distance(0), energy(0) {}

        // double return_energy(double distance);
        void set_energy(const double energy_) {this->energy = energy_;}

        void set_type(const int type_id) { this->type = type_id; }
        [[nodiscard]] int get_type() const { return type; }

        void set_UB_force_constant(const double force_constant) {this->UB_force_constant = force_constant; }
        [[nodiscard]] double get_UB_force_constant() const {return UB_force_constant; }

        void set_UB_equil_value(const double equil_value) {this->UB_force_equil = equil_value; }
        [[nodiscard]] double get_UB_equil_value() const {return UB_force_equil; }

        void set_distance_value(const double distance_) {this->distance = distance_;}

        [[nodiscard]] int get_atomA_index() const {return atomA_index; }
        [[nodiscard]] int get_atomB_index() const {return atomB_index; }
        [[nodiscard]] double return_energy(double distance) const;

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