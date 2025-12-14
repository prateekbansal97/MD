#ifndef CMapGroup_H
#define CMapGroup_H

#include <vector>
#include <utility>
class CMapGroup
{

    public:
        
        CMapGroup(const int parameter_set, const int atomA_index, const int atomB_index, const int atomC_index, const int atomD_index, const int atomE_index) :
        parameter_set(parameter_set),
        atomA_index(atomA_index),
        atomB_index(atomB_index),
        atomC_index(atomC_index),
        atomD_index(atomD_index),
        atomE_index(atomE_index),
        phi(0.0), psi(0.0), energy(0.0) {}


        static double return_energy_linear(double phi, double psi, const std::vector<double>& grid, int resolution);
        static std::pair<double, double> return_gradient_linear(double phi, double psi, const std::vector<double>& grid, int resolution);

        [[nodiscard]] double return_energy_bicubic(double phi, double psi, int resolution, const std::vector<double>& coeffs, const std::vector<double>& grid_energies) const;
        [[nodiscard]] std::pair<double, double> return_gradient_bicubic(double phi, double psi, int resolution, const std::vector<double>& coeffs, const std::vector<double>& grid_energies) const;

        void set_energy(const double energy_) {this->energy = energy_;}

        void set_parameter_set(const int set) { this->parameter_set = set; }
        [[nodiscard]] int get_parameter_set() const { return parameter_set; }


        void set_angles(const double phi_, const double psi_) {
            this->phi = phi_;
            this->psi = psi_;
        }

        [[nodiscard]] double get_angle1() const { return phi; }
        [[nodiscard]] double get_angle2() const { return psi; }

        [[nodiscard]] int get_atomA_index() const {return atomA_index; }
        [[nodiscard]] int get_atomB_index() const {return atomB_index; }
        [[nodiscard]] int get_atomC_index() const {return atomC_index; }
        [[nodiscard]] int get_atomD_index() const {return atomD_index; }
        [[nodiscard]] int get_atomE_index() const {return atomE_index; }


    private:
        int parameter_set;
        int atomA_index;
        int atomB_index;
        int atomC_index;
        int atomD_index;
        int atomE_index;
        double phi;
        double psi;
        double energy;
};
#endif //CMapGroup_H