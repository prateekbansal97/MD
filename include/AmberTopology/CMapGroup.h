#ifndef CMapGroup_H
#define CMapGroup_H


class CMapGroup
{

    public:
        
        CMapGroup(int parameter_set, int atomA_index, int atomB_index, int atomC_index, int atomD_index, int atomE_index) : 
        atomA_index(atomA_index), 
        atomB_index(atomB_index), 
        atomC_index(atomC_index),
        atomD_index(atomD_index),
        atomE_index(atomE_index),
        parameter_set(parameter_set) {}



        // double return_energy(double dihedral);

        void set_parameter_set(int set) { this->parameter_set = set; }
        const int get_parameter_set() const { return parameter_set; }


        void set_angles(double phi, double psi) {
            this->phi = phi;
            this->psi = psi;
        }

        double get_angle1() const { return phi; }
        double get_angle2() const { return psi; }

        const int get_atomA_index() const {return atomA_index; }
        const int get_atomB_index() const {return atomB_index; }
        const int get_atomC_index() const {return atomC_index; }
        const int get_atomD_index() const {return atomD_index; }
        const int get_atomE_index() const {return atomE_index; }


    private:
        int parameter_set;
        int atomA_index;
        int atomB_index;
        int atomC_index;
        int atomD_index;
        int atomE_index;
        double phi;
        double psi;
};
#endif //CMapGroup_H