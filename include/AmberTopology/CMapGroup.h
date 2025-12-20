#ifndef CMapGroup_H
#define CMapGroup_H

#include <vector>
#include <utility>

class CMapGroup
{
    public:
        
        CMapGroup(const int parameter_set, const int atomA_index, const int atomB_index, const int atomC_index, const int atomD_index, const int atomE_index) :
        parameter_set_(parameter_set),
        atomA_index_(atomA_index),
        atomB_index_(atomB_index),
        atomC_index_(atomC_index),
        atomD_index_(atomD_index),
        atomE_index_(atomE_index),
        phi_(0.0), psi_(0.0), energy_(0.0) {}



        void set_energy(const double energy) {this->energy_ = energy;}
        void set_parameter_set(const int set) { this->parameter_set_ = set; }
        void set_angles(const double phi, const double psi) {
            this->phi_ = phi;
            this->psi_ = psi;
        }

        [[nodiscard]] double get_angle1() const { return phi_; }
        [[nodiscard]] double get_angle2() const { return psi_; }
        [[nodiscard]] int get_atomA_index() const {return atomA_index_; }
        [[nodiscard]] int get_atomB_index() const {return atomB_index_; }
        [[nodiscard]] int get_atomC_index() const {return atomC_index_; }
        [[nodiscard]] int get_atomD_index() const {return atomD_index_; }
        [[nodiscard]] int get_atomE_index() const {return atomE_index_; }
        [[nodiscard]] int get_parameter_set() const { return parameter_set_; }
        [[nodiscard]] double get_energy() const { return energy_; }



    private:
        static double return_energy_linear(double phi, double psi, const std::vector<double>& grid, int resolution);
        static std::pair<double, double> return_gradient_linear(double phi, double psi, const std::vector<double>& grid, int resolution);
        [[nodiscard]] double return_energy_bicubic(double phi, double psi, int resolution, const std::vector<double>& coeffs, const std::vector<double>& grid_energies) const;
        [[nodiscard]] std::pair<double, double> return_gradient_bicubic(double phi, double psi, int resolution, const std::vector<double>& coeffs, const std::vector<double>& grid_energies) const;

        int parameter_set_;
        int atomA_index_;
        int atomB_index_;
        int atomC_index_;
        int atomD_index_;
        int atomE_index_;
        double phi_;
        double psi_;
        double energy_;

        friend class System;
};
#endif //CMapGroup_H