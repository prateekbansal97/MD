//
// Created by Prateek Bansal on 12/20/25.
//

#include "AmberTopology/AmberTopology.h"
#include <fstream>
#include "util/IO.h"

namespace md::AmberTopology
{
    [[maybe_unused]] void AmberTopology::print_atom_details(const int max_print) const
    {
        int printed = 0;
        for (auto& atom : atom_list_) {
            if (printed >= max_print) break; // Print only first 1000 atoms for brevity
            printed++;
            std::cout << "Atom Index " << printed
                      << ", Atom Name: " << atom.get_name()
                      << ", Type: " << atom.get_type()
                      << ", Charge: " << atom.get_partial_charge()
                      << ", Atomic Number: " << atom.get_atomic_number()
                      << ", Mass: " << atom.get_mass()
                      << ", Element: " << atom.get_element()
                      << ", Atom Type Index: " << atom.get_atom_type_index()
                      << ", Number of Excluded Atoms: " << atom.get_nExcluded_Atoms()
                      << ", Residue Name: " << atom.get_residue_name()
                      << ", Residue Number: " << atom.get_residue_number()
                      << std::endl;
        }
    }

    [[maybe_unused]] void AmberTopology::print_atom_details_to_file() const
    {
        // Optional: Implement writing atom details to a file if needed
        std::ofstream outfile("/Users/Prateek/Desktop/Coding/MD/atom_details.txt");
        if (!outfile) {
            std::cerr << "Error opening output file." << std::endl;
            return;
        }
        for (auto& atom : atom_list_) {
            outfile << "\""<< atom.get_name() << "\", ";
        }
    }

    [[maybe_unused]] void AmberTopology::print_UB_bonds(const int max_print)
    {
        if (HarmonicUB_list_.empty()) {
            std::cout << "No Urey-Bradley bonds parsed." << std::endl;
            return;
        }

        size_t count = HarmonicUB_list_.size();
        if (max_print > 0 && static_cast<size_t>(max_print) < count) {
            count = static_cast<size_t>(max_print);
        }

        for (size_t i = 0; i < count; i++)
        {
            const Bonded::HarmonicUB& UB_bond = HarmonicUB_list_[i];
            const int atomA_index = UB_bond.get_atomA_index();
            const int atomB_index = UB_bond.get_atomB_index();

            check_if_valid_indices("atomA_index", atom_list_, atomA_index);
            check_if_valid_indices("atomB_index", atom_list_, atomB_index);

            Atom& atomA = atom_list_[atomA_index];
            Atom& atomB = atom_list_[atomB_index];


            std::cout << "UB[" << i << "] "
                      << "type=" << UB_bond.get_type() + 1
                      << " k=" << UB_bond.get_UB_force_constant()
                      << " r0=" << UB_bond.get_UB_equil_value()
                      << " atoms: " << atomA.get_name() << " (" << atomA.get_residue_name() << atomA.get_residue_number() << ")"
                      << " - " << atomB.get_name() << " (" << atomB.get_residue_name() << atomB.get_residue_number() << ")"
                      << std::endl;
        }
    }

    void AmberTopology::process_dihedral_force_constant(const std::string& line)
    {
        for (const std::vector<std::string> entries = util::split_line_over_empty_spaces(line); const auto & entry : entries) {
            double force_constant = std::stod(entry);
            dihedral_force_constants_.push_back(force_constant);
        }
    }

    [[maybe_unused]] void AmberTopology::print_IMP_bonds(const int max_print)
    {
        if (HarmonicIMP_list_.empty()) {
            std::cout << "No Improper torsions parsed." << std::endl;
            return;
        }

        size_t count = HarmonicIMP_list_.size();
        if (max_print > 0 && static_cast<size_t>(max_print) < count) {
            count = static_cast<size_t>(max_print);
        }

        for (size_t i = 0; i < count; i++)
        {
            const Bonded::HarmonicImproper& IMP_bond = HarmonicIMP_list_[i];
            const int atomA_index = IMP_bond.get_atomA_index();
            const int atomB_index = IMP_bond.get_atomB_index();
            const int atomC_index = IMP_bond.get_atomC_index();
            const int atomD_index = IMP_bond.get_atomD_index();

            check_if_valid_indices("atom_list_", atom_list_, atomA_index, atomB_index, atomC_index, atomD_index);

            Atom& atomA = atom_list_[atomA_index];
            Atom& atomB = atom_list_[atomB_index];
            Atom& atomC = atom_list_[atomC_index];
            Atom& atomD = atom_list_[atomC_index];

            std::cout << "IMP[" << i << "] "
                      << "type=" << IMP_bond.get_type() + 1
                      << " k=" << IMP_bond.get_IMP_force_constant()
                      << " phase=" << IMP_bond.get_IMP_phase_value()
                      << " atoms: " << atomA.get_name() << " (" << atomA.get_residue_name() << atomA.get_residue_number() << ")"
                      << " - " << atomB.get_name() << " (" << atomB.get_residue_name() << atomB.get_residue_number() << ")"
                      << " - " << atomC.get_name() << " (" << atomC.get_residue_name() << atomC.get_residue_number() << ")"
                      << " - " << atomD.get_name() << " (" << atomD.get_residue_name() << atomD.get_residue_number() << ")"
                      << std::endl;
        }
    }

    [[maybe_unused]] void AmberTopology::print_bonds_without_H(const int max_print, const int start_point)
    {
        if (HarmonicBond_list_.empty()) {
            std::cout << "No bonds parsed." << std::endl;
            return;
        }

        size_t count = HarmonicBond_list_.size();
        std::cout << "Total number of bonds: " << count << std::endl;
        if (max_print > 0 && static_cast<size_t>(max_print) < count && start_point + max_print < count) {
            count = static_cast<size_t>(max_print);
        }

        for (size_t i = start_point; i < count + start_point; i++)
        {
            check_if_valid_indices("HarmonicBond_list_", HarmonicBond_list_, i);
            const Bonded::HarmonicBond& bond = HarmonicBond_list_[i];
            const int atomA_index = bond.get_atomA_index();
            const int atomB_index = bond.get_atomB_index();
            check_if_valid_indices("atom_list_", atom_list_, atomA_index, atomB_index);

            Atom& atomA = atom_list_[atomA_index];
            Atom& atomB = atom_list_[atomB_index];


            std::cout << "Bond[" << i << "] "
                      << "type=" << bond.get_type() + 1
                      << " k=" << bond.get_Bond_force_constant()
                      << " r0=" << bond.get_Bond_equil_length()
                      << " atoms: " << atomA.get_name() << " (" << atomA.get_residue_name() << atomA.get_residue_number() << ")"
                      << " - " << atomB.get_name() << " (" << atomB.get_residue_name() << atomB.get_residue_number() << ")"
                      << std::endl;
        }
    }

    [[maybe_unused]] void AmberTopology::print_bonds_including_H(const int max_print, const int start_point)
    {
        if (HarmonicBond_list_.empty()) {
            std::cout << "No bonds parsed." << std::endl;
            return;
        }

        size_t count = HarmonicBond_list_.size();
        std::cout << "Total number of bonds: " << count << std::endl;
        if (max_print > 0 && static_cast<size_t>(max_print) < count && start_point + max_print < count) {
            count = static_cast<size_t>(max_print);
        }

        for (size_t i = start_point; i < count + start_point; i++)
        {
            check_if_valid_indices("HarmonicBond_list_", HarmonicBond_list_, i);
            const Bonded::HarmonicBond& bond = HarmonicBond_list_[i];
            const int atomA_index = bond.get_atomA_index();
            const int atomB_index = bond.get_atomB_index();
            check_if_valid_indices("atom_list_", atom_list_, atomA_index, atomB_index);

            Atom& atomA = atom_list_[atomA_index];
            Atom& atomB = atom_list_[atomB_index];

            std::cout << "Bond[" << i << "] "
                      << "type=" << bond.get_type() + 1
                      << " k=" << bond.get_Bond_force_constant()
                      << " r0=" << bond.get_Bond_equil_length()
                      << " atoms: " << atomA.get_name() << " (" << atomA.get_residue_name() << atomA.get_residue_number() << ")"
                      << " - " << atomB.get_name() << " (" << atomB.get_residue_name() << atomB.get_residue_number() << ")"
                      << std::endl;
        }
    }

    [[maybe_unused]] void AmberTopology::print_angles_without_H(const int max_print, const int start_point)
    {
        if (HarmonicAngle_list_.empty()) {
            std::cout << "No angles parsed." << std::endl;
            return;
        }

        size_t count = HarmonicAngle_list_.size();
        std::cout << "Total number of angles: " << count << std::endl;
        if (max_print > 0 && static_cast<size_t>(max_print) < count && start_point + max_print < count) {
            count = static_cast<size_t>(max_print);
        }

        for (size_t i = start_point; i < count + start_point; i++)
        {
            check_if_valid_indices("HarmonicAngle_list_", HarmonicAngle_list_, i);
            const Bonded::HarmonicAngle& angle = HarmonicAngle_list_[i];
            const int atomA_index = angle.get_atomA_index();
            const int atomB_index = angle.get_atomB_index();
            const int atomC_index = angle.get_atomC_index();
            check_if_valid_indices("atom_list_", atom_list_, atomA_index, atomB_index, atomC_index);


            Atom& atomA = atom_list_[atomA_index];
            Atom& atomB = atom_list_[atomB_index];
            Atom& atomC = atom_list_[atomC_index];
            std::cout << "Angle[" << i << "] "
                      << "type=" << angle.get_type() + 1
                      << " k=" << angle.get_Angle_force_constant()
                      << " r0=" << angle.get_Angle_equil_angle()
                      << " atoms: " << atomA.get_name() << " (" << atomA.get_residue_name() << atomA.get_residue_number() << ")"
                      << " - " << atomB.get_name() << " (" << atomB.get_residue_name() << atomB.get_residue_number() << ")"
                      << " - " << atomC.get_name() << " (" << atomC.get_residue_name() << atomC.get_residue_number() << ")"
                      << std::endl;
        }
    }

    [[maybe_unused]] void AmberTopology::print_angles_including_H(const int max_print, const int start_point)
    {
        if (HarmonicAngle_list_.empty()) {
            std::cout << "No angles parsed." << std::endl;
            return;
        }

        size_t count = HarmonicAngle_list_.size();
        std::cout << "Total number of angles: " << count << std::endl;
        if (max_print > 0 && static_cast<size_t>(max_print) < count && start_point + max_print < count) {
            count = static_cast<size_t>(max_print);
        }

        for (size_t i = start_point; i < count + start_point; i++)
        {
            check_if_valid_indices("HarmonicAngle_list_", HarmonicAngle_list_, i);

            const Bonded::HarmonicAngle& angle = HarmonicAngle_list_[i];
            const int atomA_index = angle.get_atomA_index();
            const int atomB_index = angle.get_atomB_index();
            const int atomC_index = angle.get_atomC_index();
            check_if_valid_indices("atom_list_", atom_list_, atomA_index, atomB_index, atomC_index);

            Atom& atomA = atom_list_[atomA_index];
            Atom& atomB = atom_list_[atomB_index];
            Atom& atomC = atom_list_[atomC_index];
            std::cout << "Angle[" << i << "] "
                      << "type=" << angle.get_type() + 1
                      << " k=" << angle.get_Angle_force_constant()
                      << " r0=" << angle.get_Angle_equil_angle()
                      << " atoms: " << atomA.get_name() << " (" << atomA.get_residue_name() << atomA.get_residue_number() << ")"
                      << " - " << atomB.get_name() << " (" << atomB.get_residue_name() << atomB.get_residue_number() << ")"
                      << " - " << atomC.get_name() << " (" << atomC.get_residue_name() << atomC.get_residue_number() << ")"
                      << std::endl;
        }
    }

    [[maybe_unused]] void AmberTopology::print_dihedrals_without_H(const int max_print, const int start_point)
    {
        if (CosineDihedral_list_.empty()) {
            std::cout << "No dihedrals parsed." << std::endl;
            return;
        }

        size_t count = CosineDihedral_list_.size();
        if (max_print > 0 && static_cast<size_t>(max_print) < count && start_point + max_print < count) {
            count = static_cast<size_t>(max_print);
        }

        for (size_t i = start_point; i < count + start_point; i++)
        {
            check_if_valid_indices("CosineDihedral_list_", CosineDihedral_list_, i);
            const Bonded::CosineDihedral& dihedral = CosineDihedral_list_[i];
            const int atomA_index = dihedral.get_atomA_index();
            const int atomB_index = dihedral.get_atomB_index();
            const int atomC_index = dihedral.get_atomC_index();
            const int atomD_index = dihedral.get_atomD_index();
            check_if_valid_indices("atom_list_", atom_list_, atomA_index, atomB_index, atomC_index, atomD_index);

            Atom& atomA = atom_list_[atomA_index];
            Atom& atomB = atom_list_[atomB_index];
            Atom& atomC = atom_list_[atomC_index];
            Atom& atomD = atom_list_[atomC_index];
            check_if_valid_indices("atom_list_", atom_list_, atomA_index, atomB_index, atomC_index, atomD_index);

            std::cout << "Dihedral[" << i << "] "
                      << "type=" << dihedral.get_type() + 1
                      << " k=" << dihedral.get_Dihedral_force_constant()
                      << " r0=" << dihedral.get_Dihedral_phase()
                      << " atoms: " << atomA.get_name() << " (" << atomA.get_residue_name() << atomA.get_residue_number() << ")"
                      << " - " << atomB.get_name() << " (" << atomB.get_residue_name() << atomB.get_residue_number() << ")"
                      << " - " << atomC.get_name() << " (" << atomC.get_residue_name() << atomC.get_residue_number() << ")"
                      << " - " << atomD.get_name() << " (" << atomD.get_residue_name() << atomD.get_residue_number() << ")"
                      << std::endl;
        }
    }

    [[maybe_unused]] void AmberTopology::print_dihedrals_including_H(const int max_print, const int start_point)
    {
        if (CosineDihedral_list_.empty()) {
            std::cout << "No dihedrals parsed." << std::endl;
            return;
        }

        size_t count = CosineDihedral_list_.size();

        if (max_print > 0 && static_cast<size_t>(max_print) < count && start_point + max_print < count) {
            count = static_cast<size_t>(max_print);
        }

        for (size_t i = start_point; i < count + start_point; i++)
        {
            check_if_valid_indices("CosineDihedral_list_", CosineDihedral_list_, i);
            const Bonded::CosineDihedral& dihedral = CosineDihedral_list_[i];
            const int atomA_index = dihedral.get_atomA_index();
            const int atomB_index = dihedral.get_atomB_index();
            const int atomC_index = dihedral.get_atomC_index();
            const int atomD_index = dihedral.get_atomD_index();
            check_if_valid_indices("atom_list_", atom_list_, atomA_index, atomB_index, atomC_index, atomD_index);

            Atom& atomA = atom_list_[atomA_index];
            Atom& atomB = atom_list_[atomB_index];
            Atom& atomC = atom_list_[atomC_index];
            Atom& atomD = atom_list_[atomC_index];
            std::cout << "Dihedral[" << i << "] "
                      << "type=" << dihedral.get_type() + 1
                      << " k=" << dihedral.get_Dihedral_force_constant()
                      << " r0=" << dihedral.get_Dihedral_phase()
                      << " atoms: " << atomA.get_name() << " (" << atomA.get_residue_name() << atomA.get_residue_number() << ")"
                      << " - " << atomB.get_name() << " (" << atomB.get_residue_name() << atomB.get_residue_number() << ")"
                      << " - " << atomC.get_name() << " (" << atomC.get_residue_name() << atomC.get_residue_number() << ")"
                      << " - " << atomD.get_name() << " (" << atomD.get_residue_name() << atomD.get_residue_number() << ")"
                      << std::endl;
        }
    }

}
