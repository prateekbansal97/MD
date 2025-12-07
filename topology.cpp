#include <atomic>
#include <chrono>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <map>
#include <vector>

#include "topology.h"
#include "harmonicUB.h"
#include "atom.h"
#include "io.h"
#include "element.h"
#include "HarmonicBond.h"
#include "HarmonicAngle.h"
#include "CosineDihedral.h"

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
    Charmm_Urey_Bradley_Count,
    Charmm_Urey_Bradley,
    Charmm_Urey_Bradley_Force_Constant,
    Charmm_Urey_Bradley_Equil_Value,
    Dihedral_Force_Constant,
    Dihedral_Periodicity,
    Dihedral_Phase,
    Scee_Scale_Factor,
    Scnb_Scale_Factor,
    Charmm_Num_Impropers,
    Charmm_Impropers,
    Charmm_Num_Impr_Types,
    Charmm_Improper_Force_Constants,
    Charmm_Improper_Phases,
    Lennard_Jones_Acoef,
    Lennard_Jones_Bcoef,
    Lennard_Jones_14_Acoef,
    Lennard_Jones_14_Bcoef,
    Bonds_Inc_Hydrogen,
    Bonds_without_Hydrogen,
    Angles_Inc_Hydrogen,
    Angles_without_Hydrogen,
    Dihedrals_Inc_Hydrogen,
    Dihedrals_without_Hydrogen,
    Excluded_Atoms_List,
    // add more as needed
};


std::vector<unsigned long int> Topology::get_pointers() 
{
    return pointers_;
}

int Topology::read_topology(std::string filename)
{
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
                if (header == "POINTERS") current_section = Section::Pointers;
                
                else if (header == "ATOM_NAME") current_section = Section::AtomName;
                else if (header == "CHARGE") { processed_atoms_index = 0; current_section = Section::Charge; }
                else if (header == "ATOMIC_NUMBER") { processed_atoms_index = 0; current_section = Section::Atomic_Number; }
                else if (header == "MASS") { processed_atoms_index = 0; current_section = Section::Mass; }
                else if (header == "ATOM_TYPE_INDEX") { processed_atoms_index = 0; current_section = Section::Atom_Type_Index; }
                
                else if (header == "NUMBER_EXCLUDED_ATOMS") { processed_atoms_index = 0; current_section = Section::Number_Excluded_Atoms; }
                else if (header == "NONBONDED_PARM_INDEX"){ processed_atoms_index = 0; index_processed = 0; current_section = Section::Nonbonded_Parm_Index; }
                
                else if (header == "RESIDUE_LABEL") { index_processed = 0; current_section = Section::Residue_Label; }
                else if (header == "RESIDUE_POINTER") { res_num = 0; current_section = Section::Residue_Pointer; }
                
                else if (header == "BOND_FORCE_CONSTANT") { res_num = 0; current_section = Section::Bond_Force_Constant; }
                else if (header == "BOND_EQUIL_VALUE") { res_num = 0; current_section = Section::Bond_Equil_Value; }
                
                else if (header == "ANGLE_FORCE_CONSTANT") current_section = Section::Angle_Force_Constant; 
                else if (header == "ANGLE_EQUIL_VALUE") current_section = Section::Angle_Equil_Value; 
                
                else if (header == "CHARMM_UREY_BRADLEY_COUNT") current_section = Section::Charmm_Urey_Bradley_Count;
                else if (header == "CHARMM_UREY_BRADLEY") current_section = Section::Charmm_Urey_Bradley;
                else if (header == "CHARMM_UREY_BRADLEY_FORCE_CONSTANT") current_section = Section::Charmm_Urey_Bradley_Force_Constant;
                else if (header == "CHARMM_UREY_BRADLEY_EQUIL_VALUE") current_section = Section::Charmm_Urey_Bradley_Equil_Value;
                
                else if (header == "DIHEDRAL_FORCE_CONSTANT") current_section = Section::Dihedral_Force_Constant;
                else if (header == "DIHEDRAL_PERIODICITY") current_section = Section::Dihedral_Periodicity;
                else if (header == "DIHEDRAL_PHASE") current_section = Section::Dihedral_Phase;
                
                else if (header == "SCEE_SCALE_FACTOR") current_section = Section::Scee_Scale_Factor;
                else if (header == "SCNB_SCALE_FACTOR") current_section = Section::Scnb_Scale_Factor;
                
                else if (header == "CHARMM_NUM_IMPROPERS") current_section = Section::Charmm_Num_Impropers;
                else if (header == "CHARMM_IMPROPERS") current_section = Section::Charmm_Impropers;
                else if (header == "CHARMM_NUM_IMPR_TYPES") current_section = Section::Charmm_Num_Impr_Types;
                else if (header == "CHARMM_IMPROPER_FORCE_CONSTANT") current_section = Section::Charmm_Improper_Force_Constants;
                else if (header == "CHARMM_IMPROPER_PHASE") current_section = Section::Charmm_Improper_Phases;
                
                else if (header == "LENNARD_JONES_ACOEF") current_section = Section::Lennard_Jones_Acoef;
                else if (header == "LENNARD_JONES_BCOEF") current_section = Section::Lennard_Jones_Bcoef;
                else if (header == "LENNARD_JONES_14_ACOEF") current_section = Section::Lennard_Jones_14_Acoef;
                else if (header == "LENNARD_JONES_14_BCOEF") current_section = Section::Lennard_Jones_14_Bcoef;

                else if (header == "BONDS_INC_HYDROGEN") current_section = Section::Bonds_Inc_Hydrogen;
                else if (header == "BONDS_WITHOUT_HYDROGEN") current_section = Section::Bonds_without_Hydrogen;

                else if (header == "ANGLES_INC_HYDROGEN") current_section = Section::Angles_Inc_Hydrogen;
                else if (header == "ANGLES_WITHOUT_HYDROGEN") current_section = Section::Angles_without_Hydrogen;

                else if (header == "DIHEDRALS_INC_HYDROGEN") current_section = Section::Dihedrals_Inc_Hydrogen;
                else if (header == "DIHEDRALS_WITHOUT_HYDROGEN") current_section = Section::Dihedrals_without_Hydrogen;
                
                else if (header == "EXCLUDED_ATOMS_LIST") current_section = Section::Excluded_Atoms_List;
                
                else current_section = Section::None;
            }

            else {
                std::string header = extract_header(line);
                std::cout << "Skipping unhandled header: " << header << std::endl;
                continue;
            }
        }
        else {
            if (current_section == Section::Pointers) process_pointers_section(line);
            
            else if (current_section == Section::AtomName) { assign_hyperparameters(); process_atom_names_section(line); }
            else if (current_section == Section::Charge) process_charge_section(line);
            else if (current_section == Section::Atomic_Number) process_atomic_number_section(line);
            else if (current_section == Section::Mass) process_mass_section(line);
            else if (current_section == Section::Atom_Type_Index) process_atom_type_index_section(line);
            
            else if (current_section == Section::Number_Excluded_Atoms) process_number_excluded_atoms_section(line);
            else if (current_section == Section::Nonbonded_Parm_Index) process_nonbonded_parm_index_section(line);
            
            else if (current_section == Section::Residue_Label) process_residue_label_section(line);
            else if (current_section == Section::Residue_Pointer) process_residue_pointer_section(line);
            
            else if (current_section == Section::Bond_Force_Constant) process_bond_force_constant_section(line);
            else if (current_section == Section::Bond_Equil_Value) process_bond_equil_section(line);
            
            else if (current_section == Section::Angle_Force_Constant) process_angle_force_constant_section(line);
            else if (current_section == Section::Angle_Equil_Value) process_angle_equil_section(line);
            
            else if (current_section == Section::Charmm_Urey_Bradley_Count) process_charmm_urey_bradley_count_section(line);
            else if (current_section == Section::Charmm_Urey_Bradley) process_charmm_urey_bradley(line);
            else if (current_section == Section::Charmm_Urey_Bradley_Force_Constant) process_charmm_urey_bradley_force_constant(line);
            else if (current_section == Section::Charmm_Urey_Bradley_Equil_Value) process_charmm_urey_bradley_equil_value(line);
            
            else if (current_section == Section::Dihedral_Force_Constant) process_dihedral_force_constant(line);
            else if (current_section == Section::Dihedral_Periodicity) process_dihedral_periodicity(line);
            else if (current_section == Section::Dihedral_Phase) process_dihedral_phase(line);
            
            else if (current_section == Section::Scee_Scale_Factor) process_scee_scale_factor(line);
            else if (current_section == Section::Scnb_Scale_Factor) process_scnb_scale_factor(line);
            
            else if (current_section == Section::Charmm_Num_Impropers) process_charmm_num_impropers(line);
            else if (current_section == Section::Charmm_Impropers) process_charmm_impropers(line);
            else if (current_section == Section::Charmm_Num_Impr_Types) process_charmm_num_impr_types(line);
            else if (current_section == Section::Charmm_Improper_Force_Constants) process_charmm_improper_force_constant(line);
            else if (current_section == Section::Charmm_Improper_Phases) process_charmm_improper_phase(line);

            else if (current_section == Section::Lennard_Jones_Acoef) process_Lennard_Jones_Acoef(line);
            else if (current_section == Section::Lennard_Jones_Bcoef) process_Lennard_Jones_Bcoef(line);

            else if (current_section == Section::Bonds_Inc_Hydrogen) process_bonds_including_H(line);
            else if (current_section == Section::Bonds_without_Hydrogen) process_bonds_without_H(line);

            else if (current_section == Section::Angles_Inc_Hydrogen) process_angles_including_H(line);
            else if (current_section == Section::Angles_without_Hydrogen) process_angles_without_H(line);

            else if (current_section == Section::Dihedrals_Inc_Hydrogen) process_dihedrals_including_H(line);
            else if (current_section == Section::Dihedrals_without_Hydrogen) process_dihedrals_without_H(line);

            else if (current_section == Section::Excluded_Atoms_List) process_excluded_atoms_list(line);
            else continue;
        }
    }

    process_charmm_urey_bradley_assign();
    process_charmm_urey_bradley_assign_ffparams();
    process_charmm_improper_assign();
    process_improper_bradley_assign_ffparams();
    create_bonds_including_H();
    create_bonds_without_H();
    create_angles_including_H();
    create_angles_without_H();
    create_dihedrals_including_H();
    create_dihedrals_without_H();
    create_excluded_atoms_list();
    // std::cout << "excluded_atoms_list_ length: " << excluded_atoms_list_.size() << std::endl;
    // print_atom_details(200);
    // print_UB_bonds(100);
    // print_bonds_including_H(100, 87600);
    // print_IMP_bonds(100);
    // print_atom_details_to_file();
    // print_dihedrals_without_H(100, 0);
    
    file.close();
    return 0;
}

void Topology::print_atom_details(int max_print)
{
    int printed = 0;
    for (const auto& atom : atom_list_) {
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
    for (const auto& atom : atom_list_) {
        outfile << "\""<< atom.name << "\", ";
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
            auto it = get_element_from_name(atom.name);
            atom.element = it;
            atom_list_.push_back(atom);
            continue;
        }
        else {
            Atom atom(entry);
            auto it = get_element_from_name(atom.name);
            atom.element = it;
            atom_list_.push_back(atom);
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
            if (processed_atoms_index < atom_list_.size()) {
                atom_list_[processed_atoms_index].partial_charge = charge;
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
        if (processed_atoms_index < atom_list_.size()) {
            atom_list_[processed_atoms_index].atomic_number = atomic_num;
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
    
        double mass = std::stod(entries[i]);
        if (processed_atoms_index < atom_list_.size()) {
            atom_list_[processed_atoms_index].mass = mass;
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
        if (processed_atoms_index < atom_list_.size()) {
            atom_list_[processed_atoms_index].atom_type_index = type_index;
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
    
        int num_excluded = std::stoi(entries[i]);
        if (processed_atoms_index < atom_list_.size()) {
            atom_list_[processed_atoms_index].nExcluded_Atoms = num_excluded;
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
            atom_list_[i].residue_number = res_num;
            atom_list_[i].residue_name = residue_labels_[res_num - 1];
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
            atom_list_[i].residue_number = res_num;
            atom_list_[i].residue_name = residue_labels_[res_num - 1];
        }
    }
}   

void Topology::process_bond_force_constant_section(std::string& line)
{
    std::vector<std::string> entries = split_line_over_empty_spaces(line);

    for (size_t i = 0; i < entries.size(); ++i) {
    
        double force_constant = std::stod(entries[i]);
        bond_force_constants_.push_back(force_constant);
    }
}

void Topology::process_bond_equil_section(std::string& line)
{
        std::vector<std::string> entries = split_line_over_empty_spaces(line);

    for (size_t i = 0; i < entries.size(); ++i) {
    
        double force_constant = std::stod(entries[i]);
        bond_force_equils_.push_back(force_constant);
    }
}

void Topology::process_angle_force_constant_section(std::string& line)
{
    std::vector<std::string> entries = split_line_over_empty_spaces(line);

    for (size_t i = 0; i < entries.size(); ++i) {
    
        double force_constant = std::stod(entries[i]);
        angle_force_constants_.push_back(force_constant);
    }
}

void Topology::process_angle_equil_section(std::string& line)
{
        std::vector<std::string> entries = split_line_over_empty_spaces(line);

    for (size_t i = 0; i < entries.size(); ++i) {
    
        double force_constant = std::stod(entries[i]);
        angle_force_equils_.push_back(force_constant);
    }
}

void Topology::process_charmm_urey_bradley_count_section(std::string& line)
{
    std::vector<std::string> entries = split_line_over_empty_spaces(line);
    nUreyBradley_ =  static_cast<unsigned long int>(std::stoul(entries[0]));
    nTypesUreyBradley = static_cast<unsigned long int>(std::stoul(entries[1]));
}

void Topology::process_charmm_urey_bradley(std::string& line)
{
    std::vector<std::string> entries = split_line_over_empty_spaces(line);

    for (size_t i = 0; i < entries.size(); ++i) {
        int index = std::stoi(entries[i]);
        charmm_urey_bradley_indices_.push_back(index);
    }
}

void Topology::process_charmm_urey_bradley_assign()
{
    for (int i = 0; i + 3 < static_cast<int>(charmm_urey_bradley_indices_.size()); i=i+3)
    {
        int indexA = charmm_urey_bradley_indices_[i] - 1;
        int indexB = charmm_urey_bradley_indices_[i+1] - 1;
        int index_CUB_type = charmm_urey_bradley_indices_[i+2] - 1;

        // Atom atomA = atom_list_[indexA];
        // Atom atomB = atom_list_[indexB];
        HarmonicUB UB_bond = HarmonicUB(indexA, indexB);
        UB_bond.set_type(index_CUB_type);
        HarmonicUB_list_.push_back(UB_bond);
    }
}

void Topology::process_charmm_urey_bradley_force_constant(std::string& line)
{
    std::vector<std::string> entries = split_line_over_empty_spaces(line);
    for (size_t i = 0; i < entries.size(); ++i) {
    
        double force_constant = std::stod(entries[i]);
        charmm_urey_bradley_force_constants_.push_back(force_constant);
    }
}

void Topology::process_charmm_urey_bradley_equil_value(std::string& line)
{
    std::vector<std::string> entries = split_line_over_empty_spaces(line);
    for (size_t i = 0; i < entries.size(); ++i) {
    
        double equil_value = std::stod(entries[i]);
        charmm_urey_bradley_equil_values_.push_back(equil_value);
    }
}

void Topology::process_charmm_urey_bradley_assign_ffparams()
{
    for (auto& UB_bond: HarmonicUB_list_)
    {
        int type = UB_bond.get_type();
        double force_constant = charmm_urey_bradley_force_constants_[type];
        double equil_value = charmm_urey_bradley_equil_values_[type];
        UB_bond.set_UB_force_constant(force_constant);
        UB_bond.set_UB_equil_value(equil_value);
    }
}

void Topology::print_UB_bonds(int max_print)
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
        const HarmonicUB& UB_bond = HarmonicUB_list_[i];
        const int atomA_index = UB_bond.get_atomA_index();
        const int atomB_index = UB_bond.get_atomB_index();
        Atom& atomA = atom_list_[atomA_index];
        Atom& atomB = atom_list_[atomB_index];
        std::cout << "UB[" << i << "] "
                  << "type=" << UB_bond.get_type() + 1
                  << " k=" << UB_bond.get_UB_force_constant()
                  << " r0=" << UB_bond.get_UB_equil_value()
                  << " atoms: " << atomA.name << " (" << atomA.residue_name << atomA.residue_number << ")"
                  << " - " << atomB.name << " (" << atomB.residue_name << atomB.residue_number << ")"
                  << std::endl;
    }
}

void Topology::process_dihedral_force_constant(std::string& line)
{
    std::vector<std::string> entries = split_line_over_empty_spaces(line);
    for (size_t i = 0; i < entries.size(); ++i) {
    
        double force_constant = std::stod(entries[i]);
        dihedral_force_constants_.push_back(force_constant);
    }
}

void Topology::process_dihedral_periodicity(std::string& line)
{
    std::vector<std::string> entries = split_line_over_empty_spaces(line);
    for (size_t i = 0; i < entries.size(); ++i) {
    
        double periodicity = std::stod(entries[i]);
        dihedral_periodicities_.push_back(periodicity);
    }
}

void Topology::process_dihedral_phase(std::string& line)
{
    std::vector<std::string> entries = split_line_over_empty_spaces(line);
    for (size_t i = 0; i < entries.size(); ++i) {
    
        double phase = std::stod(entries[i]);
        dihedral_phases_.push_back(phase);
    }
}

void Topology::process_scee_scale_factor(std::string& line)
{
    std::vector<std::string> entries = split_line_over_empty_spaces(line);
    for (size_t i = 0; i < entries.size(); ++i) {
    
        double scee = std::stod(entries[i]);
        scee_scale_factors_.push_back(scee);
    }
}

void Topology::process_scnb_scale_factor(std::string& line)
{
    std::vector<std::string> entries = split_line_over_empty_spaces(line);
    for (size_t i = 0; i < entries.size(); ++i) {
    
        double scnb = std::stod(entries[i]);
        scnb_scale_factors_.push_back(scnb);
    }
}

void Topology::process_charmm_num_impropers(std::string& line)
{
    std::vector<std::string> entries = split_line_over_empty_spaces(line);
    nCharmmImpropers_ = std::stoi(entries[0]);
}

void Topology::process_charmm_impropers(std::string& line)
{
    std::vector<std::string> entries = split_line_over_empty_spaces(line);
    for (size_t i = 0; i < entries.size(); ++i) {
        int index = std::stoi(entries[i]);
        charmm_improper_indices_.push_back(index);   
    }
}

void Topology::process_charmm_improper_assign()
{
    for (int i = 0; i + 4 < static_cast<int>(charmm_improper_indices_.size()); i=i+5)
    {
        int indexA = charmm_improper_indices_[i] - 1;
        int indexB = charmm_improper_indices_[i+1] - 1;
        int indexC = charmm_improper_indices_[i+2] - 1;
        int indexD = charmm_improper_indices_[i+3] - 1;
        int index_CI_type = charmm_improper_indices_[i+4] - 1;

        // Atom atomA = atom_list_[indexA];
        // Atom atomB = atom_list_[indexB];
        // Atom atomC = atom_list_[indexC];
        // Atom atomD = atom_list_[indexD];

        HarmonicImproper IMP_bond = HarmonicImproper(indexA, indexB, indexC, indexD);
        IMP_bond.set_type(index_CI_type);
        HarmonicIMP_list_.push_back(IMP_bond);
    }
}

void Topology::process_charmm_num_impr_types(std::string& line)
{
    std::vector<std::string> entries = split_line_over_empty_spaces(line);
    nCharmmImproperTypes_ = std::stoi(entries[0]);
}

void Topology::process_charmm_improper_force_constant(std::string& line)
{
    std::vector<std::string> entries = split_line_over_empty_spaces(line);
    for (size_t i = 0; i < entries.size(); ++i) {
    
        double force_constant = std::stod(entries[i]);
        charmm_improper_force_constants_.push_back(force_constant);
    }
}

void Topology::process_charmm_improper_phase(std::string& line)
{
    std::vector<std::string> entries = split_line_over_empty_spaces(line);
    for (size_t i = 0; i < entries.size(); ++i) {
    
        double phase = std::stod(entries[i]);
        charmm_improper_phases_.push_back(phase);
    }
}

void Topology::process_improper_bradley_assign_ffparams()
{
    for (HarmonicImproper& IMP_bond: HarmonicIMP_list_)
    {
        int type = IMP_bond.get_type();
        double force_constant = charmm_improper_force_constants_[type];
        double phase_value = charmm_improper_phases_[type];
        IMP_bond.set_IMP_force_constant(force_constant);
        IMP_bond.set_IMP_phase_value(phase_value);
    }
}

void Topology::print_IMP_bonds(int max_print)
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
        const HarmonicImproper& IMP_bond = HarmonicIMP_list_[i];
        const int atomA_index = IMP_bond.get_atomA_index();
        const int atomB_index = IMP_bond.get_atomB_index();
        const int atomC_index = IMP_bond.get_atomC_index();
        const int atomD_index = IMP_bond.get_atomD_index();
        Atom& atomA = atom_list_[atomA_index];
        Atom& atomB = atom_list_[atomB_index];
        Atom& atomC = atom_list_[atomC_index];
        Atom& atomD = atom_list_[atomC_index];

        std::cout << "IMP[" << i << "] "
                  << "type=" << IMP_bond.get_type() + 1
                  << " k=" << IMP_bond.get_IMP_force_constant()
                  << " phase=" << IMP_bond.get_IMP_phase_value()
                  << " atoms: " << atomA.name << " (" << atomA.residue_name << atomA.residue_number << ")"
                  << " - " << atomB.name << " (" << atomB.residue_name << atomB.residue_number << ")"
                  << " - " << atomC.name << " (" << atomC.residue_name << atomC.residue_number << ")"
                  << " - " << atomD.name << " (" << atomD.residue_name << atomD.residue_number << ")"
                  << std::endl;
    }
}

void Topology::process_Lennard_Jones_Acoef(std::string& line)
{
    std::vector<std::string> entries = split_line_over_empty_spaces(line);
    for (size_t i = 0; i < entries.size(); ++i) {
    
        double Acoef = std::stod(entries[i]);
        lennard_jones_Acoefs_.push_back(Acoef);
    }
}

void Topology::process_Lennard_Jones_Bcoef(std::string& line)
{
    std::vector<std::string> entries = split_line_over_empty_spaces(line);
    for (size_t i = 0; i < entries.size(); ++i) {
    
        double Bcoef = std::stod(entries[i]);
        lennard_jones_Bcoefs_.push_back(Bcoef);
    }
}

void Topology::process_Lennard_Jones_14_Acoef(std::string& line)
{
    std::vector<std::string> entries = split_line_over_empty_spaces(line);
    for (size_t i = 0; i < entries.size(); ++i) {
    
        double Acoef = std::stod(entries[i]);
        lennard_jones_14_Acoefs_.push_back(Acoef);
    }
}

void Topology::process_Lennard_Jones_14_Bcoef(std::string& line)
{
    std::vector<std::string> entries = split_line_over_empty_spaces(line);
    for (size_t i = 0; i < entries.size(); ++i) {
    
        double Bcoef = std::stod(entries[i]);
        lennard_jones_14_Bcoefs_.push_back(Bcoef);
    }
}

void Topology::process_bonds_including_H(std::string& line) // BONDS_INC_HYDROGEN
{
    std::vector<std::string> entries = split_line_over_empty_spaces(line);
    for (size_t i = 0; i < entries.size(); ++i) {
    
        int Acoef = std::stoi(entries[i]);
        bonds_including_h_.push_back(Acoef);
    }
}

void Topology::create_bonds_including_H()
{
    for (size_t i = 0; i + 3 < bonds_including_h_.size(); i = i + 3)
    {
        int indexA = (bonds_including_h_[i])/3;
        int indexB = (bonds_including_h_[i+1])/3;
        int force_type = bonds_including_h_[i+2] - 1;

        // std::cout << "index A: " << indexA << " index B: " << indexB << " force_type: " << force_type << std::endl;
        if (force_type < 0)
        {
            std::cout << "Invalid indices";
        }

        // Atom atomA = atom_list_[indexA];
        // Atom atomB = atom_list_[indexB];

        double force_constant = bond_force_constants_[force_type];
        double equil_length = bond_force_equils_[force_type];

        HarmonicBond Bond = HarmonicBond(force_constant, equil_length, indexA, indexB, true);
        HarmonicBond_list_.push_back(Bond);
    }
}

void Topology::print_bonds_without_H(int max_print, int start_point)
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
        const HarmonicBond& bond = HarmonicBond_list_[i];
        const int atomA_index = bond.get_atomA_index();
        const int atomB_index = bond.get_atomB_index();
        Atom& atomA = atom_list_[atomA_index];
        Atom& atomB = atom_list_[atomB_index];
        std::cout << "Bond[" << i << "] "
                  << "type=" << bond.get_type() + 1
                  << " k=" << bond.get_Bond_force_constant()
                  << " r0=" << bond.get_Bond_equil_length()
                  << " atoms: " << atomA.name << " (" << atomA.residue_name << atomA.residue_number << ")"
                  << " - " << atomB.name << " (" << atomB.residue_name << atomB.residue_number << ")"
                  << std::endl;
    }
}

void Topology::process_bonds_without_H(std::string& line) // BONDS_WITHOUT_HYDROGEN
{
    std::vector<std::string> entries = split_line_over_empty_spaces(line);
    for (size_t i = 0; i < entries.size(); ++i) {
    
        int Acoef = std::stoi(entries[i]);
        bonds_without_h_.push_back(Acoef);
    }
}

void Topology::create_bonds_without_H()
{
    for (size_t i = 0; i + 3 < bonds_without_h_.size(); i = i + 3)
    {
        int indexA = (bonds_without_h_[i])/3;
        int indexB = (bonds_without_h_[i+1])/3;
        int force_type = bonds_without_h_[i+2] - 1;

        // std::cout << "index A: " << indexA << " index B: " << indexB << " force_type: " << force_type << std::endl;
        if (force_type < 0)
        {
            std::cout << "Invalid indices";
        }

        // Atom atomA = atom_list_[indexA];
        // Atom atomB = atom_list_[indexB];
        double force_constant = bond_force_constants_[force_type];
        double equil_length = bond_force_equils_[force_type];

        HarmonicBond Bond = HarmonicBond(force_constant, equil_length, indexA, indexB, false);
        HarmonicBond_list_.push_back(Bond);
    }
}

void Topology::print_bonds_including_H(int max_print, int start_point)
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
        const HarmonicBond& bond = HarmonicBond_list_[i];
        const int atomA_index = bond.get_atomA_index();
        const int atomB_index = bond.get_atomB_index();
        Atom& atomA = atom_list_[atomA_index];
        Atom& atomB = atom_list_[atomB_index];
        std::cout << "Bond[" << i << "] "
                  << "type=" << bond.get_type() + 1
                  << " k=" << bond.get_Bond_force_constant()
                  << " r0=" << bond.get_Bond_equil_length()
                  << " atoms: " << atomA.name << " (" << atomA.residue_name << atomA.residue_number << ")"
                  << " - " << atomB.name << " (" << atomB.residue_name << atomB.residue_number << ")"
                  << std::endl;
    }
}

void Topology::process_angles_including_H(std::string& line) // BONDS_INC_HYDROGEN
{
    std::vector<std::string> entries = split_line_over_empty_spaces(line);
    for (size_t i = 0; i < entries.size(); ++i) {
    
        int Acoef = std::stoi(entries[i]);
        angles_including_h_.push_back(Acoef);
    }
}

void Topology::create_angles_including_H()
{
    for (size_t i = 0; i + 4 < angles_including_h_.size(); i = i + 4)
    {
        int indexA = (angles_including_h_[i])/3;
        int indexB = (angles_including_h_[i+1])/3;
        int indexC = (angles_including_h_[i+2])/3;
        int force_type = angles_including_h_[i+3] - 1;

        // std::cout << "index A: " << indexA << " index B: " << indexB << " force_type: " << force_type << std::endl;
        if (force_type < 0)
        {
            std::cout << "Invalid indices";
        }

        // Atom atomA = atom_list_[indexA];
        // Atom atomB = atom_list_[indexB];
        // Atom atomC = atom_list_[indexC];
        double force_constant = angle_force_constants_[force_type];
        double equil_angle = angle_force_equils_[force_type];

        HarmonicAngle Angle = HarmonicAngle(force_constant, equil_angle, indexA, indexB, indexC, true);
        HarmonicAngle_list_.push_back(Angle);
    }
}

void Topology::print_angles_without_H(int max_print, int start_point)
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
        const HarmonicAngle& angle = HarmonicAngle_list_[i];
        const int atomA_index = angle.get_atomA_index();
        const int atomB_index = angle.get_atomB_index();
        const int atomC_index = angle.get_atomC_index();
        Atom& atomA = atom_list_[atomA_index];
        Atom& atomB = atom_list_[atomB_index];
        Atom& atomC = atom_list_[atomC_index];
        std::cout << "Angle[" << i << "] "
                  << "type=" << angle.get_type() + 1
                  << " k=" << angle.get_Angle_force_constant()
                  << " r0=" << angle.get_Angle_equil_angle()
                  << " atoms: " << atomA.name << " (" << atomA.residue_name << atomA.residue_number << ")"
                  << " - " << atomB.name << " (" << atomB.residue_name << atomB.residue_number << ")"
                  << " - " << atomC.name << " (" << atomC.residue_name << atomC.residue_number << ")"
                  << std::endl;
    }
}

void Topology::process_angles_without_H(std::string& line) // BONDS_WITHOUT_HYDROGEN
{
    std::vector<std::string> entries = split_line_over_empty_spaces(line);
    for (size_t i = 0; i < entries.size(); ++i) {
    
        int Acoef = std::stoi(entries[i]);
        angles_without_h_.push_back(Acoef);
    }
}

void Topology::create_angles_without_H()
{
    for (size_t i = 0; i + 4 < angles_without_h_.size(); i = i + 4)
    {
        int indexA = (angles_without_h_[i])/3;
        int indexB = (angles_without_h_[i+1])/3;
        int indexC = (angles_without_h_[i+2])/3;
        int force_type = angles_without_h_[i+3] - 1;

        // std::cout << "index A: " << indexA << " index B: " << indexB << " force_type: " << force_type << std::endl;
        if (force_type < 0)
        {
            std::cout << "Invalid indices";
        }

        // Atom atomA = atom_list_[indexA];
        // Atom atomB = atom_list_[indexB];
        // Atom atomC = atom_list_[indexC];
        double force_constant = angle_force_constants_[force_type];
        double equil_angle = angle_force_equils_[force_type];

        HarmonicAngle Angle = HarmonicAngle(force_constant, equil_angle, indexA, indexB, indexC, false);
        HarmonicAngle_list_.push_back(Angle);
    }
}

void Topology::print_angles_including_H(int max_print, int start_point)
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
        const HarmonicAngle& angle = HarmonicAngle_list_[i];
        const int atomA_index = angle.get_atomA_index();
        const int atomB_index = angle.get_atomB_index();
        const int atomC_index = angle.get_atomC_index();
        Atom& atomA = atom_list_[atomA_index];
        Atom& atomB = atom_list_[atomB_index];
        Atom& atomC = atom_list_[atomC_index];
        std::cout << "Angle[" << i << "] "
                  << "type=" << angle.get_type() + 1
                  << " k=" << angle.get_Angle_force_constant()
                  << " r0=" << angle.get_Angle_equil_angle()
                  << " atoms: " << atomA.name << " (" << atomA.residue_name << atomA.residue_number << ")"
                  << " - " << atomB.name << " (" << atomB.residue_name << atomB.residue_number << ")"
                  << " - " << atomC.name << " (" << atomC.residue_name << atomC.residue_number << ")"
                  << std::endl;
    }
}

void Topology::process_dihedrals_including_H(std::string& line) // DIHEDRALS_INC_HYDROGEN
{
    std::vector<std::string> entries = split_line_over_empty_spaces(line, false);
    for (size_t i = 0; i < entries.size(); ++i) {
    
        int Acoef = std::stoi(entries[i]);
        dihedrals_including_h_.push_back(Acoef);
    }
}

void Topology::create_dihedrals_including_H()
{
    for (size_t i = 0; i + 5 < dihedrals_including_h_.size(); i = i + 5)
    {
        
        bool exclude_14 = false;
        bool improper = false;
        int indexA = (dihedrals_including_h_[i])/3;
        int indexB = (dihedrals_including_h_[i+1])/3;
        int indexC = (dihedrals_including_h_[i+2])/3;
        
        if (indexC < 0)
        {
            exclude_14 = true;
            indexC = std::abs(indexC);
        } 
        int indexD = (dihedrals_including_h_[i+3])/3;
        if (indexD < 0)
        {
            improper = true;
            indexD = std::abs(indexD);
        }

        int total_atoms = atom_list_.size();

        if (indexA > total_atoms || indexB > total_atoms  || indexC > total_atoms || indexD > total_atoms )
        {
            std::cerr << "Invalid indices, indices greater than total number of atoms";
        }

        if (indexA < 0 || indexB < 0  || indexC < 0 || indexD < 0 )
        {
            std::cerr << "Invalid indices, less than 0";
        }
        int force_type = dihedrals_including_h_[i+4] - 1;
        
        if (force_type < 0 || force_type > dihedral_force_constants_.size() - 1)
        {
            std::cout << "Invalid indices";
        }
        // Atom atomA = atom_list_[indexA];
        // Atom atomB = atom_list_[indexB];
        // Atom atomC = atom_list_[indexC];
        // Atom atomD = atom_list_[indexD];
        
        double force_constant = dihedral_force_constants_[force_type];
        double phase = dihedral_phases_[force_type];
        double periodicity = dihedral_periodicities_[force_type];
        
        CosineDihedral Dihedral = CosineDihedral(force_constant, phase, periodicity, indexA, indexB, indexC, indexD, true, improper, exclude_14);
        CosineDihedral_list_.push_back(Dihedral);
    }
}

void Topology::print_dihedrals_without_H(int max_print, int start_point)
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
        const CosineDihedral& dihedral = CosineDihedral_list_[i];
        const int atomA_index = dihedral.get_atomA_index();
        const int atomB_index = dihedral.get_atomB_index();
        const int atomC_index = dihedral.get_atomC_index();
        const int atomD_index = dihedral.get_atomD_index();
        Atom& atomA = atom_list_[atomA_index];
        Atom& atomB = atom_list_[atomB_index];
        Atom& atomC = atom_list_[atomC_index];
        Atom& atomD = atom_list_[atomC_index];
        std::cout << "Dihedral[" << i << "] "
                  << "type=" << dihedral.get_type() + 1
                  << " k=" << dihedral.get_Dihedral_force_constant()
                  << " r0=" << dihedral.get_Dihedral_phase()
                  << " atoms: " << atomA.name << " (" << atomA.residue_name << atomA.residue_number << ")"
                  << " - " << atomB.name << " (" << atomB.residue_name << atomB.residue_number << ")"
                  << " - " << atomC.name << " (" << atomC.residue_name << atomC.residue_number << ")"
                  << " - " << atomD.name << " (" << atomD.residue_name << atomD.residue_number << ")"
                  << std::endl;
    }
}

void Topology::process_dihedrals_without_H(std::string& line) // BONDS_WITHOUT_HYDROGEN
{
    std::vector<std::string> entries = split_line_over_empty_spaces(line);
    for (size_t i = 0; i < entries.size(); ++i) {
    
        int Acoef = std::stoi(entries[i]);
        dihedrals_without_h_.push_back(Acoef);
    }
}

void Topology::create_dihedrals_without_H()
{
    for (size_t i = 0; i + 5 < dihedrals_without_h_.size(); i = i + 5)
    {
        bool exclude_14 = false;
        bool improper = false;
        int indexA = (dihedrals_without_h_[i])/3;
        int indexB = (dihedrals_without_h_[i+1])/3;
        int indexC = (dihedrals_without_h_[i+2])/3;
        if (indexC < 0)
        {
            exclude_14 = true;
            indexC = std::abs(indexC);
        } 
        int indexD = (dihedrals_without_h_[i+3])/3;
        if (indexD < 0)
        {
            improper = true;
            indexD = std::abs(indexD);
        }
        int force_type = dihedrals_without_h_[i+4] - 1;

        if (force_type < 0)
        {
            std::cout << "Invalid indices";
        }

        // Atom atomA = atom_list_[indexA];
        // Atom atomB = atom_list_[indexB];
        // Atom atomC = atom_list_[indexC];
        // Atom atomD = atom_list_[indexD];

        double force_constant = dihedral_force_constants_[force_type];
        double phase = dihedral_phases_[force_type];
        double periodicity = dihedral_periodicities_[force_type];

        CosineDihedral Dihedral = CosineDihedral(force_constant, phase, periodicity, indexA, indexB, indexC, indexD, false, improper, exclude_14);
        CosineDihedral_list_.push_back(Dihedral);
    }
}

void Topology::print_dihedrals_including_H(int max_print, int start_point)
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
        const CosineDihedral& dihedral = CosineDihedral_list_[i];
        const int atomA_index = dihedral.get_atomA_index();
        const int atomB_index = dihedral.get_atomB_index();
        const int atomC_index = dihedral.get_atomC_index();
        const int atomD_index = dihedral.get_atomD_index();
        Atom& atomA = atom_list_[atomA_index];
        Atom& atomB = atom_list_[atomB_index];
        Atom& atomC = atom_list_[atomC_index];
        Atom& atomD = atom_list_[atomC_index];
        std::cout << "Dihedral[" << i << "] "
                  << "type=" << dihedral.get_type() + 1
                  << " k=" << dihedral.get_Dihedral_force_constant()
                  << " r0=" << dihedral.get_Dihedral_phase()
                  << " atoms: " << atomA.name << " (" << atomA.residue_name << atomA.residue_number << ")"
                  << " - " << atomB.name << " (" << atomB.residue_name << atomB.residue_number << ")"
                  << " - " << atomC.name << " (" << atomC.residue_name << atomC.residue_number << ")"
                  << " - " << atomD.name << " (" << atomD.residue_name << atomD.residue_number << ")"
                  << std::endl;
    }
}

void Topology::process_excluded_atoms_list(std::string& line)
{
    std::vector<std::string> entries = split_line_over_empty_spaces(line);
    for (size_t i = 0; i < entries.size(); ++i) {
        int excluded = std::stoi(entries[i]);
        excluded_atoms_list_.push_back(excluded);
    }
}

void Topology::create_excluded_atoms_list()
{
    size_t start = 0;
    
    // Safety check: Ensure the list is large enough
    // This prevents a crash if the file is truncated
    size_t total_expected = excluded_atoms_list_.size(); 

    for (auto& atom: atom_list_)
    {
        unsigned long int excluded = atom.nExcluded_Atoms;

        // Bounds check
        if (start + excluded > total_expected) {
            std::cerr << "Error: Excluded atoms list is shorter than expected!" << std::endl;
            break; 
        }

        for (size_t p = start; p < start + excluded; p++)
        {
            int amber_index = excluded_atoms_list_[p];

            // AMBER index 0 is a placeholder/null. 
            // If we find 0, we shouldn't add -1 to the list.
            if (amber_index > 0) {
                atom.excluded_atoms.push_back(amber_index - 1);
            }
        }
        start += excluded;
    }
}

