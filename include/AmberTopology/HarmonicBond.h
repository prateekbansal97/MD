#ifndef HARMONICBOND_H
#define HARMONICBOND_H




class HarmonicBond
{

    public:
        HarmonicBond(const double Bond_force_constant, const double Bond_Equil) :
            Bond_force_constant(Bond_force_constant),
            Bond_Equil(Bond_Equil), atomA_index(-1),
            atomB_index(-1),
            type(0),
            isH(false), distance(0), energy(0) {}

        HarmonicBond(const double Bond_force_constant, const double Bond_Equil, const int atomA_index, const int atomB_index, const bool isH) :
            Bond_force_constant(Bond_force_constant),
            Bond_Equil(Bond_Equil),
            atomA_index(atomA_index),
            atomB_index(atomB_index), type(0),
            isH(isH), distance(0), energy(0) {}

        HarmonicBond(const int atomA_index, const int atomB_index) :
            Bond_force_constant(0.0),
            Bond_Equil(0.0),
            atomA_index(atomA_index),
            atomB_index(atomB_index), type(0),
            isH(false), distance(0), energy(0) {}

        [[nodiscard]] double return_energy(double distance) const;
        void set_energy(const double energy_) { this->energy = energy_;}

        void set_type(const int type_id) { this->type = type_id; }
        [[nodiscard]] int get_type() const { return type; }

        void set_Bond_force_constant(const double force_constant) {this->Bond_force_constant = force_constant; }
        [[nodiscard]] double get_Bond_force_constant() const {return Bond_force_constant; }

        void set_Bond_equil_length(const double equil_length) {this->Bond_Equil = equil_length; }
        [[nodiscard]] double get_Bond_equil_length() const {return Bond_Equil; }

        void set_distance(const double distance_) {this->distance = distance_;}

        [[nodiscard]] int get_atomA_index() const {return atomA_index; }
        [[nodiscard]] int get_atomB_index() const {return atomB_index; }

        [[nodiscard]] bool get_isH() const {return isH;}


    private:
        double Bond_force_constant;
        double Bond_Equil;
        int atomA_index;
        int atomB_index;
        int type;
        bool isH;
        double distance;
        double energy;
};
#endif //HARMONICBOND_H