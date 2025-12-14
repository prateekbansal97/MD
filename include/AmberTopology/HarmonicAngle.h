#ifndef HARMONICANGLE_H
#define HARMONICANGLE_H

class HarmonicAngle
{

    public:
        HarmonicAngle(const double Angle_force_constant, const double Angle_Equil) :
            Angle_force_constant(Angle_force_constant),
            Angle_Equil(Angle_Equil), atomA_index(-1),
            atomB_index(-1),
            atomC_index(-1),
            type(0),
            isH(false), angle(0), energy(0) {}

        HarmonicAngle(const double Angle_force_constant, const double Angle_Equil, const int atomA_index, const int atomB_index, const int atomC_index, const bool isH) :
            Angle_force_constant(Angle_force_constant),
            Angle_Equil(Angle_Equil),
            atomA_index(atomA_index),
            atomB_index(atomB_index), atomC_index(atomC_index),
            type(0), isH(isH), angle(0), energy(0) {}

        HarmonicAngle(const int atomA_index, const int atomB_index, const int atomC_index) :
            Angle_force_constant(0.0),
            Angle_Equil(0.0),
            atomA_index(atomA_index),
            atomB_index(atomB_index), atomC_index(atomC_index), type(0), isH(false),
            angle(0), energy(0) {}

        [[nodiscard]] double return_energy(double angle) const;
        void set_energy(const double energy_) {this->energy = energy_;}

        void set_type(const int type_id) { this->type = type_id; }
        [[nodiscard]] int get_type() const { return type; }

        void set_Angle_force_constant(const double force_constant) {this->Angle_force_constant = force_constant; }
        [[nodiscard]] double get_Angle_force_constant() const {return Angle_force_constant; }

        void set_Angle_equil_angle(const double equil_angle) {this->Angle_Equil = equil_angle; }
        [[nodiscard]] double get_Angle_equil_angle() const {return Angle_Equil; }

        void set_angle(const double angle_) { this->angle = angle_;}

        [[nodiscard]] int get_atomA_index() const {return atomA_index; }
        [[nodiscard]] int get_atomB_index() const {return atomB_index; }
        [[nodiscard]] int get_atomC_index() const {return atomC_index; }

        [[nodiscard]] bool get_isH() const {return isH;}


    private:
        double Angle_force_constant;
        double Angle_Equil;
        int atomA_index;
        int atomB_index;
        int atomC_index;
        int type;
        bool isH;
        double angle;
        double energy;

};
#endif //HARMONICANGLE_H