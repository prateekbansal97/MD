#include <chrono>
#include <iostream>
#include <fstream>
#include <string>
#include "topology.h"
#include <vector>
#include "atom.h"
#include "io.h"
#include <map>
#include "element.h"

enum class Section {
    None,
    Pointers,
    AtomName,
    Charge,
    Atomic_Number,
    Mass,
    Atom_Type_Index,
    Number_Excluded_Atoms,
    Nonbonded_Parm_Index,
    Residue_Label,
    Residue_Pointer,
    Bond_Force_Constant,
    Bond_Equil_Value,
    Angle_Force_Constant,
    Angle_Equil_Value,
    // add more as needed
};


std::vector<unsigned long int> Topology::get_pointers() 
{
    return pointers_;
}

int Topology::read_topology()
{
    std::string filename = "/Users/Prateek/Desktop/Coding/MD/CB1_apo_assym.prmtop";
    std::ifstream file(filename);

    if (!file) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return 1;
    }


    std::string line;
    Section current_section = Section::None;

    while (std::getline(file, line)) {


        if (check_if_line_empty(line)) {
            continue;
        }


        if (check_first_char(line, '%'))
        {
            if (check_if_line_starts_with_string(line, "%FORMAT"))
            {
                continue;
            }

            else if (check_if_line_starts_with_string(line, "%COMMENT"))
            {
                std::cout << "[COMMENT] " << extract_header(line) << std::endl;
                continue;
            }


            else if (check_if_line_starts_with_string(line, "%FLAG"))
            {
                std::string header = extract_header(line);
                if (header == "POINTERS")
                {
                    current_section = Section::Pointers;
                }
                else if (header == "ATOM_NAME")
                {
                    current_section = Section::AtomName;
                }
                else if (header == "CHARGE")
                {
                    processed_atoms_index = 0;
                    current_section = Section::Charge;
                }
                else if (header == "ATOMIC_NUMBER")
                {
                    processed_atoms_index = 0;
                    current_section = Section::Atomic_Number;
                }
                else if (header == "MASS")
                {
                    processed_atoms_index = 0;
                    current_section = Section::Mass;
                }
                else if (header == "ATOM_TYPE_INDEX")
                {
                    processed_atoms_index = 0;
                    current_section = Section::Atom_Type_Index;
                }
                else if (header == "NUMBER_EXCLUDED_ATOMS")
                {
                    processed_atoms_index = 0;
                    current_section = Section::Number_Excluded_Atoms;
                }
                else if (header == "NONBONDED_PARM_INDEX")
                {
                    processed_atoms_index = 0;
                    index_processed = 0;
                    current_section = Section::Nonbonded_Parm_Index;
                }
                else if (header == "RESIDUE_LABEL")
                {
                    index_processed = 0;
                    current_section = Section::Residue_Label;
                }
                else if (header == "RESIDUE_POINTER")
                {
                    res_num = 0;
                    current_section = Section::Residue_Pointer;
                }
                else if (header == "BOND_FORCE_CONSTANT")
                {
                    res_num = 0;
                    current_section = Section::Bond_Force_Constant;
                }
                else if (header == "BOND_EQUIL_VALUE")
                {
                    res_num = 0;
                    current_section = Section::Bond_Equil_Value;
                }
                else if (header == "ANGLE_FORCE_CONSTANT")
                {
                    current_section = Section::Angle_Force_Constant;
                }
                else if (header == "ANGLE_EQUIL_VALUE")
                {
                    current_section = Section::Angle_Equil_Value;
                }
                else
                {   
                    current_section = Section::None;
                }

            }

            else {
                std::string header = extract_header(line);
                std::cout << "Skipping unhandled header: " << header << std::endl;
                continue;
            }
        }
        else {
            if (current_section == Section::Pointers) {
                process_pointers_section(line);
            }
            else if (current_section == Section::AtomName) {
                assign_hyperparameters();
                process_atom_names_section(line);
            }
            else if (current_section == Section::Charge) {
                process_charge_section(line);
            }
            else if (current_section == Section::Atomic_Number) {
                process_atomic_number_section(line);
            }
            else if (current_section == Section::Mass) {
                process_mass_section(line);
            }
            else if (current_section == Section::Atom_Type_Index) {
                process_atom_type_index_section(line);
            }
            else if (current_section == Section::Number_Excluded_Atoms) {
                process_number_excluded_atoms_section(line);
            }
            else if (current_section == Section::Nonbonded_Parm_Index) {
                process_nonbonded_parm_index_section(line);
            }
            else if (current_section == Section::Residue_Label) {
                process_residue_label_section(line);
            }
            else if (current_section == Section::Residue_Pointer) {
                process_residue_pointer_section(line); 
            }
            else if (current_section == Section::Bond_Force_Constant)
            {
                process_bond_force_constant_section(line);
            }
            else if (current_section == Section::Bond_Equil_Value)
            {
                process_bond_equil_section(line);
            }
            else if (current_section == Section::Angle_Force_Constant)
            {
                process_angle_force_constant_section(line);
            }
            else if (current_section == Section::Angle_Equil_Value)
            {
                process_angle_equil_section(line);
            }
            else {
                // Skip unhandled sections
                continue;
            }
        }
    }


    print_atom_details(200);
    // print_atom_details_to_file();

    file.close();
    return 0;
}

void Topology::print_atom_details(int max_print)
{
    int printed = 0;
    for (const auto& atom : atom_list) {
        if (printed >= max_print) break; // Print only first 1000 atoms for brevity
        printed++;
        std::cout << "Atom Index " << printed
                  << ", Atom Name: " << atom.name
                  << ", Type: " << atom.type
                  << ", Charge: " << atom.partial_charge
                  << ", Atomic Number: " << atom.atomic_number
                  << ", Mass: " << atom.mass
                  << ", Element: " << atom.element
                  << ", Atom Type Index: " << atom.atom_type_index
                  << ", Number of Excluded Atoms: " << atom.nExcluded_Atoms
                  << ", Residue Name: " << atom.residue_name
                  << ", Residue Number: " << atom.residue_number
                  << std::endl;
    }
}

void Topology::print_atom_details_to_file()
{
    // Optional: Implement writing atom details to a file if needed
    std::ofstream outfile("/Users/Prateek/Desktop/Coding/MD/atom_details.txt");
    if (!outfile) {
        std::cerr << "Error opening output file." << std::endl;
        return;
    }
    for (const auto& atom : atom_list) {
        outfile << "\""<< atom.name << "\", ";
                // << ", Type: " << atom.type
                // << ", Charge: " << atom.partial_charge
                // << ", Atomic Number: " << atom.atomic_number
                // << ", Mass: " << atom.mass
                // << ", Element: " << atom.element
                // << ", Atom Type Index: " << atom.atom_type_index
                // << ", Number of Excluded Atoms: " << atom.nExcluded_Atoms
                // << ", Nonbonded Parm Index: " << atom.nonbonded_parm_index
                // << std::endl;
    }
}


void Topology::process_pointers_section(std::string& line)
{
    std::vector<std::string> entries = split_line_over_empty_spaces(line);
    for (const std::string& entry : entries) {
        try {
            unsigned long int value = static_cast<unsigned long int>(std::stoi(entry));
            pointers_.push_back(value);
        } catch (const std::invalid_argument& e) {
            std::cerr << "Invalid number format: " << entry << std::endl;
        } catch (const std::out_of_range& e) {
            std::cerr << "Number out of range: " << entry << std::endl;
        }
    }
}

void Topology::assign_hyperparameters()
{
    if (pointers_.size() < 29) {
        std::cerr << "Insufficient POINTERS data: expected at least 29 entries, got " << pointers_.size() << std::endl;
        return;
    }
    nAtoms_ = pointers_[0];
    nTypes_ = pointers_[1];
    nBondHs_ = pointers_[2];
    nBondAs_ = pointers_[3];
    nAnglesH_ = pointers_[4];
    nAnglesA_ = pointers_[5];
    nTorsionsH_ = pointers_[6];
    nTorsionsA_ = pointers_[7];
    nexclusions_ = pointers_[10];
    nresidues_ = pointers_[11];
    nconstraint_bonds = pointers_[12];
    nconstraint_angles = pointers_[13];
    nconstraint_torsions_ = pointers_[14];
    nunique_bond_types_ = pointers_[15];
    nunique_angle_types_ = pointers_[16];
    nunique_torsion_types_ = pointers_[17];
    nsolty_terms_ = pointers_[18];
    nphb_types_ = pointers_[19];
    ifpert_ = pointers_[20];
    nperturbed_bonds_ = pointers_[21];
    nperturbed_angles_ = pointers_[22];
    nperturbed_torsions_ = pointers_[23];
    mperturbed_bonds_ = pointers_[24];
    mperturbed_angles_ = pointers_[25];
    mperturbed_torsions_ = pointers_[26];
    ifbox_ = pointers_[27];
    nmaxresidue_atoms_ = pointers_[28];

    // Initialize nbmatrix_ to nTypes_ x nTypes_
    if (nTypes_ > 0) {
        try {
            nbmatrix_.assign(static_cast<size_t>(nTypes_), std::vector<unsigned long int>(static_cast<size_t>(nTypes_), 0));
        } catch (const std::exception& e) {
            std::cerr << "Failed to allocate nbmatrix_ with size " << nTypes_ << "x" << nTypes_ << ": " << e.what() << std::endl;
        }
    } else {
        std::cerr << "Warning: nTypes_ is zero; nbmatrix_ not initialized." << std::endl;
    }
}

void Topology::process_atom_names_section(std::string& line)
{
    std::vector<std::string> entries = split_line_fixed_length(line, 4);
    for (const auto& entry : entries) {
        if (check_if_whitespace_in_string(entry)) {
            Atom atom(remove_whitespaces_from_string(entry));
            // Use the map in element.h to assign element
            auto it = element_map.find(atom.name);
            if (it != element_map.end()) {
                atom.element = it->second;
            } else {
                atom.element = atom.name.substr(0, 1); // Default to first character
            }
            atom_list.push_back(atom);
            continue;
        }
        else {
            Atom atom(entry);
            auto it = element_map.find(atom.name);
            if (it != element_map.end()) {
                atom.element = it->second;
            } else {
                atom.element = atom.name.substr(0, 1); // Default to first character
            }
            atom_list.push_back(atom);
        }
    }
}

void Topology::process_charge_section(std::string& line)
{
    std::vector<std::string> entries = split_line_over_empty_spaces(line);
    for (size_t i = 0; i < entries.size(); ++i) {
        try {
            long double charge = std::stold(entries[i]);

            // divide by sqrt(332.0716D0) to convert to elementary charge units
            charge /= 18.2206899283247L;

            // std::cout << "Assigning charge " << charge << " to atom index " << i << std::endl;
            if (processed_atoms_index < atom_list.size()) {
                atom_list[processed_atoms_index].partial_charge = charge;
                processed_atoms_index++;
            } else {
                std::cerr << "Charge index out of range: " << i << std::endl;
            }
        } catch (const std::invalid_argument& e) {
            std::cerr << "Invalid charge format: " << entries[i] << std::endl;
        } catch (const std::out_of_range& e) {
            std::cerr << "Charge out of range: " << entries[i] << std::endl;
        }
    }
}

void Topology::process_atomic_number_section(std::string& line)
{
    std::vector<std::string> entries = split_line_over_empty_spaces(line);
    for (size_t i = 0; i < entries.size(); ++i) {
    
        int atomic_num = std::stoi(entries[i]);
        if (processed_atoms_index < atom_list.size()) {
            atom_list[processed_atoms_index].atomic_number = atomic_num;
            processed_atoms_index++;
        } else {
            std::cerr << "Atomic number index out of range: " << i << std::endl;        
    }
}
}

void Topology::process_mass_section(std::string& line)
{
    std::vector<std::string> entries = split_line_over_empty_spaces(line);  
    for (size_t i = 0; i < entries.size(); ++i) {
    
        float mass = std::stof(entries[i]);
        if (processed_atoms_index < atom_list.size()) {
            atom_list[processed_atoms_index].mass = mass;
            processed_atoms_index++;
        } else {
            std::cerr << "Mass index out of range: " << i << std::endl;        
    }
}
}

void Topology::process_atom_type_index_section(std::string& line)
{
    std::vector<std::string> entries = split_line_over_empty_spaces(line);  
    for (size_t i = 0; i < entries.size(); ++i) {
    
        unsigned long int type_index = static_cast<unsigned long int>(std::stoul(entries[i]));
        if (processed_atoms_index < atom_list.size()) {
            atom_list[processed_atoms_index].atom_type_index = type_index;
            processed_atoms_index++;
        } else {
            std::cerr << "Atom type index out of range: " << i << std::endl;        
    }
}
}

void Topology::process_number_excluded_atoms_section(std::string& line)
{
    std::vector<std::string> entries = split_line_over_empty_spaces(line);  
    for (size_t i = 0; i < entries.size(); ++i) {
    
        unsigned long int num_excluded = static_cast<unsigned long int>(std::stoul(entries[i]));
        if (processed_atoms_index < atom_list.size()) {
            atom_list[processed_atoms_index].nExcluded_Atoms = num_excluded;
            processed_atoms_index++;
        } else {
            std::cerr << "Number of excluded atoms index out of range: " << i << std::endl;        
    }
}
}

void Topology::process_nonbonded_parm_index_section(std::string& line)
{
    std::vector<std::string> entries = split_line_over_empty_spaces(line);
    for (size_t i = 0; i < entries.size(); ++i) {        
        unsigned long int nb_parm_index = static_cast<unsigned long int>(std::stoul(entries[i]));

        if (nTypes_ == 0) {
            std::cerr << "Nonbonded parm index section encountered but nTypes_ is zero. Skipping." << std::endl;
            return;
        }

        int row = static_cast<int>(index_processed / nTypes_);
        int col = static_cast<int>(index_processed % nTypes_);
        
        
        if (row < 0 || col < 0 || static_cast<size_t>(row) >= nbmatrix_.size() || static_cast<size_t>(col) >= nbmatrix_[row].size()) {
            std::cerr << "Index out of bounds for nbmatrix_: (" << row << ", " << col << ") - nbmatrix_ size (" << nbmatrix_.size() << ", " << (nbmatrix_.empty() ? 0 : nbmatrix_[0].size()) << ")" << std::endl;
        } else {
            nbmatrix_[row][col] = nb_parm_index;
        }
        index_processed++;
    }
}

void Topology::process_residue_label_section(std::string& line)
{
    std::vector<std::string> entries = split_line_fixed_length(line, 4);
    for (const auto& entry : entries) {
        residue_labels_.push_back(remove_whitespaces_from_string(entry));
    }
}

void Topology::process_residue_pointer_section(std::string& line)
{
    // Residue pointers are values that indicate the starting atom index for each residue.
    std::vector<std::string> entries = split_line_over_empty_spaces(line);

    if (res_num > 0)
    {
        res_num += 1;
        // int res1 = last_residue_num;
        unsigned long int res2 = static_cast<unsigned long int>(std::stoul(entries[0]));
        for (unsigned long int i = last_residue_num - 1; i < res2 - 1; i++)
        {
            atom_list[i].residue_number = res_num;
            atom_list[i].residue_name = residue_labels_[res_num - 1];
        }
    }

    for (size_t i = 0; i < entries.size() - 1; ++i)
    {

        unsigned long int res1 = static_cast<unsigned long int>(std::stoul(entries[i]));
        unsigned long int res2 = static_cast<unsigned long int>(std::stoul(entries[i+1]));

        if (i == entries.size() - 2)
        {
            last_residue_num = res2;
        }

        res_num += 1;

        for (unsigned long int i = res1 - 1; i < res2 - 1; i++)
        {
            atom_list[i].residue_number = res_num;
            atom_list[i].residue_name = residue_labels_[res_num - 1];
        }
    }
}   

void Topology::process_bond_force_constant_section(std::string& line)
{
    std::vector<std::string> entries = split_line_over_empty_spaces(line);

    for (size_t i = 0; i < entries.size(); ++i) {
    
        float force_constant = std::stof(entries[i]);
        bond_force_constants_.push_back(force_constant);
}
}

void Topology::process_bond_equil_section(std::string& line)
{
        std::vector<std::string> entries = split_line_over_empty_spaces(line);

    for (size_t i = 0; i < entries.size(); ++i) {
    
        float force_constant = std::stof(entries[i]);
        bond_force_equils_.push_back(force_constant);
}
}

void Topology::process_angle_force_constant_section(std::string& line)
{
    std::vector<std::string> entries = split_line_over_empty_spaces(line);

    for (size_t i = 0; i < entries.size(); ++i) {
    
        float force_constant = std::stof(entries[i]);
        angle_force_constants_.push_back(force_constant);
}
}

void Topology::process_angle_equil_section(std::string& line)
{
        std::vector<std::string> entries = split_line_over_empty_spaces(line);

    for (size_t i = 0; i < entries.size(); ++i) {
    
        float force_constant = std::stof(entries[i]);
        angle_force_equils_.push_back(force_constant);
}
}






