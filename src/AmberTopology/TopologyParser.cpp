//
// Created by Prateek Bansal on 12/21/25.
//

#include "AmberTopology/AmberTopology.h"
#include "util/IO.h"
#include "ChemicalEntity/Element.h"

namespace md::AmberTopology
{
    void AmberTopology::process_pointers_section(const std::string& line)
    {
        for (const std::vector<std::string> entries = util::split_line_over_empty_spaces(line); const std::string& entry : entries) {
            try {
                auto value = static_cast<unsigned long int>(std::stoi(entry));
                pointers_.push_back(value);
            } catch (const std::invalid_argument& e) {
                std::cerr << "Invalid number format: " << entry << " (" << e.what() << ")\n";
            } catch (const std::out_of_range& e) {
                std::cerr << "Number out of range: " << entry << " (" << e.what() << ")\n";
            }
        }
    }

    void AmberTopology::process_atom_names_section(const std::string& line)
    {
        for (const std::vector<std::string> entries = util::split_line_fixed_length(line, 4); const auto& entry : entries) {
            if (util::check_if_whitespace_in_string(entry)) {
                Atom atom(util::remove_whitespaces_from_string(entry));
                auto it = get_element_from_name(atom.get_name());
                atom.set_element(it);
                atom_list_.push_back(atom);
                continue;
            }
            else {
                Atom atom(entry);
                auto it = get_element_from_name(atom.get_name());
                atom.set_element(it);
                atom_list_.push_back(atom);
            }
        }
    }

    void AmberTopology::process_charge_section(const std::string& line)
    {
        const std::vector<std::string> entries = util::split_line_over_empty_spaces(line);
        for (size_t i = 0; i < entries.size(); ++i) {
            try {
                long double charge = std::stold(entries[i]);

                // divide by sqrt(332.0716D0) to convert to elementary charge units
                charge /= 18.2206899283247L;

                if (processed_atoms_index_ < atom_list_.size()) {
                    atom_list_[processed_atoms_index_].set_partial_charge(charge);
                    processed_atoms_index_++;
                } else {
                    std::cerr << "Charge index out of range: " << i << std::endl;
                }
            } catch (const std::invalid_argument& e) {
                std::cerr << "Invalid charge format: " << entries[i] << " (" << e.what() << ")\n";
            } catch (const std::out_of_range& e) {
                std::cerr << "Charge out of range: " << entries[i] << " (" << e.what() << ")\n";
            }
        }
    }

    void AmberTopology::process_atomic_number_section(const std::string& line)
    {
        const std::vector<std::string> entries = util::split_line_over_empty_spaces(line);
        for (size_t i = 0; i < entries.size(); ++i) {

            const int atomic_num = std::stoi(entries[i]);
            if (processed_atoms_index_ < atom_list_.size()) {
                atom_list_[processed_atoms_index_].set_atomic_number(atomic_num);
                processed_atoms_index_++;
            } else {
                std::cerr << "Atomic number index out of range: " << i << std::endl;
            }
        }
    }

    void AmberTopology::process_mass_section(const std::string& line)
    {
        const std::vector<std::string> entries = util::split_line_over_empty_spaces(line);
        for (size_t i = 0; i < entries.size(); ++i) {

            const double mass = std::stod(entries[i]);
            if (processed_atoms_index_ < atom_list_.size()) {
                atom_list_[processed_atoms_index_].set_mass(mass);
                processed_atoms_index_++;
            } else {
                std::cerr << "Mass index out of range: " << i << std::endl;
            }
        }
    }

    void AmberTopology::process_atom_type_index_section(const std::string& line)
    {
        const std::vector<std::string> entries = util::split_line_over_empty_spaces(line);
        for (size_t i = 0; i < entries.size(); ++i) {

            const auto type_index = static_cast<unsigned long int>(std::stoul(entries[i]));
            if (processed_atoms_index_ < atom_list_.size()) {
                atom_list_[processed_atoms_index_].set_atom_type_index(type_index);
                processed_atoms_index_++;
            } else {
                std::cerr << "Atom type index out of range: " << i << std::endl;
            }
        }
    }

    void AmberTopology::process_number_excluded_atoms_section(const std::string& line)
    {
        const std::vector<std::string> entries = util::split_line_over_empty_spaces(line);
        for (size_t i = 0; i < entries.size(); ++i) {

            const int num_excluded = std::stoi(entries[i]);
            if (processed_atoms_index_ < atom_list_.size()) {
                atom_list_[processed_atoms_index_].set_nExcluded_Atoms(num_excluded);
                processed_atoms_index_++;
            } else {
                std::cerr << "Number of excluded atoms index out of range: " << i << std::endl;
            }
        }
    }

    void AmberTopology::process_nonbonded_parm_index_section(const std::string& line)
    {
        for (const std::vector<std::string> entries = util::split_line_over_empty_spaces(line); const auto & entry : entries) {
            const auto nb_parm_index = static_cast<unsigned long int>(std::stoul(entry));

            if (nTypes_ == 0) {
                std::cerr << "Nonbonded parm index section encountered but nTypes_ is zero. Skipping." << std::endl;
                return;
            }

            const int row = static_cast<int>(index_processed_ / nTypes_);


            if (const int col = static_cast<int>(index_processed_ % nTypes_); row < 0 || col < 0 || static_cast<size_t>(row) >= nbmatrix_.size() || static_cast<size_t>(col) >= nbmatrix_[row].size()) {
                std::cerr << "Index out of bounds for nbmatrix_: (" << row << ", " << col << ") - nbmatrix_ size (" << nbmatrix_.size() << ", " << (nbmatrix_.empty() ? 0 : nbmatrix_[0].size()) << ")" << std::endl;
            } else {
                nbmatrix_[row][col] = nb_parm_index;
            }
            index_processed_++;
        }
    }

    void AmberTopology::process_residue_label_section(const std::string& line)
    {
        for (const std::vector<std::string> entries = util::split_line_fixed_length(line, 4); const auto& entry : entries) {
            residue_labels_.push_back(util::remove_whitespaces_from_string(entry));
        }
    }

    void AmberTopology::process_residue_pointer_section(const std::string& line)
    {
        // Residue pointers are values that indicate the starting atom index for each residue.
        const std::vector<std::string> entries = util::split_line_over_empty_spaces(line);

        if (res_num_ > 0)
        {
            res_num_ += 1;
            const auto res2 = static_cast<unsigned long int>(std::stoul(entries[0]));
            for (unsigned long int i = last_residue_num_ - 1; i < res2 - 1; i++)
            {
                atom_list_[i].set_residue_number(res_num_);
                atom_list_[i].set_residue_name(residue_labels_[res_num_ - 1]);
            }
        }

        for (size_t i = 0; i < entries.size() - 1; ++i)
        {

            const auto res1 = static_cast<unsigned long int>(std::stoul(entries[i]));
            const auto res2 = static_cast<unsigned long int>(std::stoul(entries[i+1]));

            if (i == entries.size() - 2)
            {
                last_residue_num_ = res2;
            }

            res_num_ += 1;

            for (unsigned long int j = res1 - 1; j < res2 - 1; j++)
            {
                atom_list_[j].set_residue_number(res_num_);
                atom_list_[j].set_residue_name(residue_labels_[res_num_ - 1]);
            }
        }
    }

    void AmberTopology::process_bond_force_constant_section(const std::string& line)
    {
        for (const std::vector<std::string> entries = util::split_line_over_empty_spaces(line); const auto & entry : entries) {

            double force_constant = std::stod(entry);
            bond_force_constants_.push_back(force_constant);
        }
    }

    void AmberTopology::process_bond_equil_section(const std::string& line)
    {
        for (const std::vector<std::string> entries = util::split_line_over_empty_spaces(line); const auto & entry : entries) {

            double force_constant = std::stod(entry);
            bond_force_equils_.push_back(force_constant);
        }
    }

    void AmberTopology::process_angle_force_constant_section(const std::string& line)
    {
        for (const std::vector<std::string> entries = util::split_line_over_empty_spaces(line); const auto & entry : entries) {

            double force_constant = std::stod(entry);
            angle_force_constants_.push_back(force_constant);
        }
    }

    void AmberTopology::process_angle_equil_section(const std::string& line)
    {
        for (const std::vector<std::string> entries = util::split_line_over_empty_spaces(line); const auto & entry : entries) {

            double force_constant = std::stod(entry);
            angle_force_equils_.push_back(force_constant);
        }
    }

    void AmberTopology::process_charmm_urey_bradley_count_section(const std::string& line)
    {
        const std::vector<std::string> entries = util::split_line_over_empty_spaces(line);
        nUreyBradley_ =  static_cast<unsigned long int>(std::stoul(entries[0]));
        nTypesUreyBradley_ = static_cast<unsigned long int>(std::stoul(entries[1]));
    }

    void AmberTopology::process_charmm_urey_bradley(const std::string& line)
    {
        for (const std::vector<std::string> entries = util::split_line_over_empty_spaces(line); const auto & entry : entries) {
            int index = std::stoi(entry);
            charmm_urey_bradley_indices_.push_back(index);
        }
    }

    void AmberTopology::process_charmm_urey_bradley_force_constant(const std::string& line)
    {
        for (const std::vector<std::string> entries = util::split_line_over_empty_spaces(line); const auto & entry : entries) {

            double force_constant = std::stod(entry);
            charmm_urey_bradley_force_constants_.push_back(force_constant);
        }
    }

    void AmberTopology::process_charmm_urey_bradley_equil_value(const std::string& line)
    {
        for (const std::vector<std::string> entries = util::split_line_over_empty_spaces(line); const auto & entry : entries) {

            double equil_value = std::stod(entry);
            charmm_urey_bradley_equil_values_.push_back(equil_value);
        }
    }

    void AmberTopology::process_dihedral_periodicity(const std::string& line)
    {
        for (const std::vector<std::string> entries = util::split_line_over_empty_spaces(line); const auto & entry : entries) {

            double periodicity = std::stod(entry);
            dihedral_periodicities_.push_back(periodicity);
        }
    }

    void AmberTopology::process_dihedral_phase(const std::string& line)
    {
        for (const std::vector<std::string> entries = util::split_line_over_empty_spaces(line); const auto & entry : entries) {

            double phase = std::stod(entry);
            dihedral_phases_.push_back(phase);
        }
    }

    void AmberTopology::process_scee_scale_factor(const std::string& line)
    {
        for (const std::vector<std::string> entries = util::split_line_over_empty_spaces(line); const auto & entry : entries) {

            double scee = std::stod(entry);
            scee_scale_factors_.push_back(scee);
        }
    }

    void AmberTopology::process_scnb_scale_factor(const std::string& line)
    {
        for (const std::vector<std::string> entries = util::split_line_over_empty_spaces(line); const auto & entry : entries) {

            double scnb = std::stod(entry);
            scnb_scale_factors_.push_back(scnb);
        }
    }

    void AmberTopology::process_charmm_num_impropers(const std::string& line)
    {
        const std::vector<std::string> entries = util::split_line_over_empty_spaces(line);
        nCharmmImpropers_ = std::stoi(entries[0]);
    }

    void AmberTopology::process_charmm_impropers(const std::string& line)
    {
        for (const std::vector<std::string> entries = util::split_line_over_empty_spaces(line); const auto & entry : entries) {
            int index = std::stoi(entry);
            charmm_improper_indices_.push_back(index);
        }
    }

    void AmberTopology::process_charmm_num_impr_types(const std::string& line)
    {
        const std::vector<std::string> entries = util::split_line_over_empty_spaces(line);
        nCharmmImproperTypes_ = std::stoi(entries[0]);
    }

    void AmberTopology::process_charmm_improper_force_constant(const std::string& line)
    {
        for (const std::vector<std::string> entries = util::split_line_over_empty_spaces(line); const auto & entry : entries) {

            double force_constant = std::stod(entry);
            charmm_improper_force_constants_.push_back(force_constant);
        }
    }

    void AmberTopology::process_charmm_improper_phase(const std::string& line)
    {
        for (const std::vector<std::string> entries = util::split_line_over_empty_spaces(line); const auto & entry : entries) {

            double phase = std::stod(entry);
            charmm_improper_phases_.push_back(phase);
        }
    }

    void AmberTopology::process_Lennard_Jones_Acoef(const std::string& line)
    {
        for (const std::vector<std::string> entries = util::split_line_over_empty_spaces(line); const auto & entry : entries) {

            double Acoef = std::stod(entry);
            lennard_jones_Acoefs_.push_back(Acoef);
        }
    }

    void AmberTopology::process_Lennard_Jones_Bcoef(const std::string& line)
    {
        for (const std::vector<std::string> entries = util::split_line_over_empty_spaces(line); const auto & entry : entries) {

            double Bcoef = std::stod(entry);
            lennard_jones_Bcoefs_.push_back(Bcoef);
        }
    }

    void AmberTopology::process_Lennard_Jones_14_Acoef(const std::string& line)
    {
        for (const std::vector<std::string> entries = util::split_line_over_empty_spaces(line); const auto & entry : entries) {

            double Acoef = std::stod(entry);
            lennard_jones_14_Acoefs_.push_back(Acoef);
        }
    }

    void AmberTopology::process_Lennard_Jones_14_Bcoef(const std::string& line)
    {
        for (const std::vector<std::string> entries = util::split_line_over_empty_spaces(line); const auto & entry : entries) {

            double Bcoef = std::stod(entry);
            lennard_jones_14_Bcoefs_.push_back(Bcoef);
        }
    }

    void AmberTopology::process_bonds_including_H(const std::string& line) // BONDS_INC_HYDROGEN
    {
        for (const std::vector<std::string> entries = util::split_line_over_empty_spaces(line); const auto & entry : entries) {

            int Acoef = std::stoi(entry);
            bonds_including_h_.push_back(Acoef);
        }
    }

    void AmberTopology::process_bonds_without_H(const std::string& line) // BONDS_WITHOUT_HYDROGEN
    {
        for (const std::vector<std::string> entries = util::split_line_over_empty_spaces(line); const auto & entry : entries) {

            int Acoef = std::stoi(entry);
            bonds_without_h_.push_back(Acoef);
        }
    }

    void AmberTopology::process_angles_including_H(const std::string& line) // BONDS_INC_HYDROGEN
    {
        for (const std::vector<std::string> entries = util::split_line_over_empty_spaces(line); const auto & entry : entries) {

            int Acoef = std::stoi(entry);
            angles_including_h_.push_back(Acoef);
        }
    }

    void AmberTopology::process_angles_without_H(const std::string& line) // BONDS_WITHOUT_HYDROGEN
    {
        for (const std::vector<std::string> entries = util::split_line_over_empty_spaces(line); const auto & entry : entries) {

            int Acoef = std::stoi(entry);
            angles_without_h_.push_back(Acoef);
        }
    }


    void AmberTopology::process_dihedrals_including_H(const std::string& line) // DIHEDRALS_INC_HYDROGEN
    {
        for (const std::vector<std::string> entries = util::split_line_over_empty_spaces(line, false); const auto & entry : entries) {

            int Acoef = std::stoi(entry);
            dihedrals_including_h_.push_back(Acoef);
        }
    }

    void AmberTopology::process_dihedrals_without_H(const std::string& line) // BONDS_WITHOUT_HYDROGEN
    {
        for (const std::vector<std::string> entries = util::split_line_over_empty_spaces(line); const auto & entry : entries) {

            int Acoef = std::stoi(entry);
            dihedrals_without_h_.push_back(Acoef);
        }
    }

    void AmberTopology::process_excluded_atoms_list(const std::string& line)
    {
        for (const std::vector<std::string> entries = util::split_line_over_empty_spaces(line); const auto & entry : entries) {
            int excluded = std::stoi(entry);
            excluded_atoms_list_.push_back(excluded);
        }
    }

    void AmberTopology::process_amber_atom_type(const std::string& line)
    {
        for (const std::vector<std::string> entries = util::split_line_fixed_length(line, 4); const auto& entry : entries) {
            amber_atom_type_.push_back(util::remove_whitespaces_from_string(entry));
        }
    }

    void AmberTopology::process_Charmm_Cmap_Count(const std::string& line)
    {
        const std::vector<std::string> entries = util::split_line_over_empty_spaces(line);
        //%COMMENT Number of CMAP terms, number of unique CMAP parameters
        nCmap_ = std::stoul(entries[0]);
        nCmap_unique_ = std::stoul(entries[1]);
    }

    void AmberTopology::process_Charmm_Cmap_Resolution(const std::string& line)
    //=CHARMM_CMAP_RESOLUTION
    {
        for (const std::vector<std::string> entries = util::split_line_over_empty_spaces(line); const auto & entry : entries) {
            int index = std::stoi(entry);
            charmm_cmap_resolution_.push_back(index);
        }
    }

    void AmberTopology::process_Charmm_Cmap_parameter_01(const std::string& line)
    //CHARMM_CMAP_PARAMETER_01
    {
        for (const std::vector<std::string> entries = util::split_line_over_empty_spaces(line); const auto & entry : entries) {
            double index = std::stod(entry);
            charmm_cmap_parameter_01_.push_back(index);
        }
    }

    void AmberTopology::process_Charmm_Cmap_parameter_02(const std::string& line)
    //CHARMM_CMAP_PARAMETER_02
    {
        for (const std::vector<std::string> entries = util::split_line_over_empty_spaces(line); const auto & entry : entries) {
            double index = std::stod(entry);
            charmm_cmap_parameter_02_.push_back(index);
        }
    }

    void AmberTopology::process_Charmm_Cmap_parameter_03(const std::string& line)
    //CHARMM_CMAP_PARAMETER_03
    {
        for (const std::vector<std::string> entries = util::split_line_over_empty_spaces(line); const auto & entry : entries) {
            double index = std::stod(entry);
            charmm_cmap_parameter_03_.push_back(index);
        }
    }

    void AmberTopology::process_Charmm_Cmap_parameter_04(const std::string& line)
    //CHARMM_CMAP_PARAMETER_04
    {
        for (const std::vector<std::string> entries = util::split_line_over_empty_spaces(line); const auto & entry : entries) {
            double index = std::stod(entry);
            charmm_cmap_parameter_04_.push_back(index);
        }
    }

    void AmberTopology::process_Charmm_Cmap_parameter_05(const std::string& line)
    //CHARMM_CMAP_PARAMETER_05
    {
        for (const std::vector<std::string> entries = util::split_line_over_empty_spaces(line); const auto & entry : entries) {
            double index = std::stod(entry);
            charmm_cmap_parameter_05_.push_back(index);
        }
    }

    void AmberTopology::process_Charmm_Cmap_Index(const std::string& line)
    {
        for (const std::vector<std::string> entries = util::split_line_over_empty_spaces(line); const auto & entry : entries) {
            int index = std::stoi(entry);
            charmm_cmap_index_.push_back(index);
        }
    }

    void AmberTopology::process_solvent_pointers(const std::string& line)
    {
        const std::vector<std::string> entries = util::split_line_over_empty_spaces(line);
        protein_res_ = std::stoul(entries[0]);
        nsolvent_mols_ = std::stoul(entries[1]);
        natoms_per_solvent_ = std::stoul(entries[2]);
    }

    void AmberTopology::process_atoms_per_molecule(const std::string& line)
    {
        for (const std::vector<std::string> entries = util::split_line_over_empty_spaces(line); const auto & entry : entries) {
            int index = std::stoi(entry);
            atoms_per_molecule_.push_back(index);
        }
    }

    void AmberTopology::process_box_dimensions(const std::string& line)
    {
        const std::vector<std::string> entries = util::split_line_over_empty_spaces(line);
        box_beta_ = std::stod(entries[0]);
        box_x_ = std::stod(entries[1]);
        box_y_ = std::stod(entries[2]);
        box_z_ = std::stod(entries[3]);
    }

    void AmberTopology::process_radius_set(const std::string& line)
    {
        radius_set_ = line;
    }

    void AmberTopology::process_radii(const std::string& line)
    {
        for (const std::vector<std::string> entries = util::split_line_over_empty_spaces(line); const auto & entry : entries) {
            double radius = std::stod(entry);
            radii_.push_back(radius);
        }
    }

    void AmberTopology::process_screen(const std::string& line)
    {
        for (const std::vector<std::string> entries = util::split_line_over_empty_spaces(line); const auto & entry : entries) {
            double screenr = std::stod(entry);
            screen_.push_back(screenr);
        }
    }

    void AmberTopology::process_ipol(const std::string& line)
    {
        const std::vector<std::string> entries = util::split_line_over_empty_spaces(line);
        polarizable_ = std::stoi(entries[0]);
    }
}