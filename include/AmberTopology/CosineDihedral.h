#ifndef COSINEDIHEDRAL_H
#define COSINEDIHEDRAL_H


// #include "atom.h"


class CosineDihedral
{

    public:
        CosineDihedral(const double Dihedral_force_constant, const double Dihedral_Phase, const double Dihedral_Periodicity,
                        const bool improper = false, const bool exclude_14 = false) :
        Dihedral_force_constant(Dihedral_force_constant),
        Dihedral_Phase(Dihedral_Phase),
        Dihedral_Periodicity(Dihedral_Periodicity), 
        atomA_index(-1),
        atomB_index(-1),
        atomC_index(-1),
        atomD_index(-1),
        type(0),
        isH(false),
        improper(improper),
        exclude_14(exclude_14) {}
        
        CosineDihedral(const double Dihedral_force_constant, const double Dihedral_Phase, const double Dihedral_Periodicity,
                        const int atomA_index, const int atomB_index, const int atomC_index, const int atomD_index,
                        const bool isH, const bool improper = false, const bool exclude_14 = false) :
        Dihedral_force_constant(Dihedral_force_constant),
        Dihedral_Phase(Dihedral_Phase),
        Dihedral_Periodicity(Dihedral_Periodicity),
        atomA_index(atomA_index), 
        atomB_index(atomB_index), 
        atomC_index(atomC_index),
        atomD_index(atomD_index),
        type(0), 
        isH(isH),
        improper(improper),
        exclude_14(exclude_14) {}

        CosineDihedral(const int atomA_index, const int atomB_index, const int atomC_index, const int atomD_index, const bool improper = false, const bool exclude_14 = false) :
        Dihedral_force_constant(0.0),
        Dihedral_Phase(0.0),
        Dihedral_Periodicity(0.0), 
        atomA_index(atomA_index), 
        atomB_index(atomB_index), 
        atomC_index(atomC_index),
        atomD_index(atomD_index),
        type(0), 
        isH(false),
        improper(improper),
        exclude_14(exclude_14) {}

        [[nodiscard]] double return_energy(double dihedral) const;
        void set_energy(const double energ) {this->energy = energ;}

        void set_type(const int type_id) { this->type = type_id; }
        [[nodiscard]] int get_type() const { return type; }

        void set_Dihedral_force_constant(const double force_constant) {this->Dihedral_force_constant = force_constant; }
        [[nodiscard]] double get_Dihedral_force_constant() const {return Dihedral_force_constant; }

        void set_Dihedral_phase(const double phase) {this->Dihedral_Phase = phase; }
        [[nodiscard]] double get_Dihedral_phase() const {return Dihedral_Phase; }

        void set_Dihedral_Periodicity(const double Periodicity) {this->Dihedral_Periodicity = Periodicity; }
        [[nodiscard]] double get_Dihedral_Periodicity() const {return Dihedral_Periodicity; }

        void set_cosine_dihedral(const double dihedra) {this->dihedral = dihedra;}

        [[nodiscard]] int get_atomA_index() const {return atomA_index; }
        [[nodiscard]] int get_atomB_index() const {return atomB_index; }
        [[nodiscard]] int get_atomC_index() const {return atomC_index; }
        [[nodiscard]] int get_atomD_index() const {return atomD_index; }
        

        [[nodiscard]] bool get_isH() const {return isH;}

        [[nodiscard]] bool get_improper() const {return improper;}
        [[nodiscard]] bool get_exclude_14() const {return exclude_14;}
    
    private:
        double Dihedral_force_constant;
        double Dihedral_Phase;
        double Dihedral_Periodicity;
        int atomA_index;
        int atomB_index;
        int atomC_index;
        int atomD_index;
        int type;
        bool isH;
        bool improper;
        bool exclude_14;
        double dihedral{};
        double energy{};
};
#endif //COSINEDIHEDRAL_H