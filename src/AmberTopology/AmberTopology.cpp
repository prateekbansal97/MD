#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <vector>
// #include <algorithm>

#include "AmberTopology/AmberTopology.h"
// #include "Bonded/harmonicUB.h"
#include "ChemicalEntity/Atom.h"
#include "util/IO.h"
// #include "AmberTopology/Element.h"
// #include "Bonded/HarmonicBond.h"
// #include "Bonded/HarmonicAngle.h"
// #include "Bonded/CosineDihedral.h"
// #include "Bonded/CMapGroup.h"
#include "ChemicalEntity/Molecule.h"

namespace md::AmberTopology
{
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
        Amber_Atom_Type,
        Charmm_Cmap_Count,
        Charmm_Cmap_Resolution,
        Charmm_Cmap_Parameter_01,
        Charmm_Cmap_Parameter_02,
        Charmm_Cmap_Parameter_03,
        Charmm_Cmap_Parameter_04,
        Charmm_Cmap_Parameter_05,
        Charmm_Cmap_Index,
        Solvent_Pointers,
        Atoms_per_Molecule,
        Box_Dimensions,
        Radius_Set,
        Radii,
        Screen,
        Ipol
    };


    [[maybe_unused]] std::vector<unsigned long int> AmberTopology::get_pointers()
    {
        return pointers_;
    }

    AmberTopology AmberTopology::read_topology_coordinates(const std::string& parmtop_path, const std::string& coords_path)
    {
        AmberTopology topology;
        std::ifstream parmtop = util::open_file(parmtop_path);
        topology.read_topology(parmtop);
        if (!coords_path.empty()) {
            std::ifstream coords = util::open_file(coords_path);
            topology.read_coords(coords);
        }
        else
        {
            std::cout << "Coordinates not supplied. Coordinates are presumed to be 0. \n";
            topology.coordinates_.assign(topology.atom_list_.size()*3, 0.0);
        }
        return topology;
    }

    int AmberTopology::read_topology(std::ifstream& parmtop)
    {
        std::string line;
        std::getline(parmtop, line);
        if (!util::check_if_line_starts_with_string(line, "%VERSION"))
        {
            std::cerr << "Invalid File Format. File does not start with %VERSION.";
            return 1;
        }

        Section current_section = Section::None;
        while (std::getline(parmtop, line)) {


            if (util::check_if_line_empty(line)) {
                continue;
            }


            if (util::check_first_char(line, '%'))
            {
                if (util::check_if_line_starts_with_string(line, "%FORMAT") || util::check_if_line_starts_with_string(line, "%COMMENT"))  continue;
                else if (util::check_if_line_starts_with_string(line, "%FLAG"))
                {
                    if (std::string header = util::extract_header(line); header == "POINTERS") current_section = Section::Pointers;

                    else if (header == "ATOM_NAME") current_section = Section::AtomName;
                    else if (header == "CHARGE") { processed_atoms_index_ = 0; current_section = Section::Charge; }
                    else if (header == "ATOMIC_NUMBER") { processed_atoms_index_ = 0; current_section = Section::Atomic_Number; }
                    else if (header == "MASS") { processed_atoms_index_ = 0; current_section = Section::Mass; }
                    else if (header == "ATOM_TYPE_INDEX") { processed_atoms_index_ = 0; current_section = Section::Atom_Type_Index; }

                    else if (header == "NUMBER_EXCLUDED_ATOMS") { processed_atoms_index_ = 0; current_section = Section::Number_Excluded_Atoms; }
                    else if (header == "NONBONDED_PARM_INDEX"){ processed_atoms_index_ = 0; index_processed_ = 0; current_section = Section::Nonbonded_Parm_Index; }

                    else if (header == "RESIDUE_LABEL") { index_processed_ = 0; current_section = Section::Residue_Label; }
                    else if (header == "RESIDUE_POINTER") { res_num_ = 0; current_section = Section::Residue_Pointer; }

                    else if (header == "BOND_FORCE_CONSTANT") { res_num_ = 0; current_section = Section::Bond_Force_Constant; }
                    else if (header == "BOND_EQUIL_VALUE") { res_num_ = 0; current_section = Section::Bond_Equil_Value; }

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
                    else if (header == "AMBER_ATOM_TYPE") current_section = Section::Amber_Atom_Type;

                    else if ((header == "CHARMM_CMAP_COUNT") || (header == "CMAP_COUNT")) current_section = Section::Charmm_Cmap_Count;
                    else if ((header == "CHARMM_CMAP_RESOLUTION") || (header == "CMAP_RESOLUTION")) current_section = Section::Charmm_Cmap_Resolution;
                    else if ((header == "CHARMM_CMAP_PARAMETER_01") || (header == "CMAP_PARAMETER_01")) current_section = Section::Charmm_Cmap_Parameter_01;
                    else if ((header == "CHARMM_CMAP_PARAMETER_02") || (header == "CMAP_PARAMETER_02")) current_section = Section::Charmm_Cmap_Parameter_02;
                    else if ((header == "CHARMM_CMAP_PARAMETER_03") || (header == "CMAP_PARAMETER_03")) current_section = Section::Charmm_Cmap_Parameter_03;
                    else if ((header == "CHARMM_CMAP_PARAMETER_04") || (header == "CMAP_PARAMETER_04")) current_section = Section::Charmm_Cmap_Parameter_04;
                    else if ((header == "CHARMM_CMAP_PARAMETER_05") || (header == "CMAP_PARAMETER_05")) current_section = Section::Charmm_Cmap_Parameter_05;
                    else if ((header == "CHARMM_CMAP_INDEX") || (header == "CMAP_INDEX")) current_section = Section::Charmm_Cmap_Index;

                    else if (header == "SOLVENT_POINTERS") current_section = Section::Solvent_Pointers;

                    else if (header == "ATOMS_PER_MOLECULE") current_section = Section::Atoms_per_Molecule;

                    else if (header == "BOX_DIMENSIONS") current_section = Section::Box_Dimensions;
                    else if (header == "RADIUS_SET") current_section = Section::Radius_Set;
                    else if (header == "RADII") current_section = Section::Radii;
                    else if (header == "SCREEN") current_section = Section::Screen;
                    else if (header == "IPOL") current_section = Section::Ipol;

                    else current_section = Section::None;
                }

                else {
                    std::string header = util::extract_header(line);
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
                else if (current_section == Section::Lennard_Jones_14_Acoef) process_Lennard_Jones_14_Acoef(line);
                else if (current_section == Section::Lennard_Jones_14_Bcoef) process_Lennard_Jones_14_Bcoef(line);

                else if (current_section == Section::Bonds_Inc_Hydrogen) process_bonds_including_H(line);
                else if (current_section == Section::Bonds_without_Hydrogen) process_bonds_without_H(line);

                else if (current_section == Section::Angles_Inc_Hydrogen) process_angles_including_H(line);
                else if (current_section == Section::Angles_without_Hydrogen) process_angles_without_H(line);

                else if (current_section == Section::Dihedrals_Inc_Hydrogen) process_dihedrals_including_H(line);
                else if (current_section == Section::Dihedrals_without_Hydrogen) process_dihedrals_without_H(line);

                else if (current_section == Section::Excluded_Atoms_List) process_excluded_atoms_list(line);
                else if (current_section == Section::Amber_Atom_Type) process_amber_atom_type(line);

                else if (current_section == Section::Charmm_Cmap_Count) process_Charmm_Cmap_Count(line);
                else if (current_section == Section::Charmm_Cmap_Resolution) process_Charmm_Cmap_Resolution(line);
                else if (current_section == Section::Charmm_Cmap_Parameter_01) process_Charmm_Cmap_parameter_01(line);
                else if (current_section == Section::Charmm_Cmap_Parameter_02) process_Charmm_Cmap_parameter_02(line);
                else if (current_section == Section::Charmm_Cmap_Parameter_03) process_Charmm_Cmap_parameter_03(line);
                else if (current_section == Section::Charmm_Cmap_Parameter_04) process_Charmm_Cmap_parameter_04(line);
                else if (current_section == Section::Charmm_Cmap_Parameter_05) process_Charmm_Cmap_parameter_05(line);
                else if (current_section == Section::Charmm_Cmap_Index) process_Charmm_Cmap_Index(line);

                else if (current_section == Section::Solvent_Pointers) process_solvent_pointers(line);

                else if (current_section == Section::Atoms_per_Molecule) process_atoms_per_molecule(line);

                else if (current_section == Section::Box_Dimensions) process_box_dimensions(line);
                else if (current_section == Section::Radius_Set) process_radius_set(line);
                else if (current_section == Section::Radii) process_radii(line);
                else if (current_section == Section::Screen) process_screen(line);
                else if (current_section == Section::Ipol) process_ipol(line);

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
        create_Charmm_Cmap_Index_assign();
        create_Cmap_Coefficient_Matrix_bicubic_spline();
        build_lj14_pairlist();
        create_molecules();

        // print_atom_details(200);
        // print_UB_bonds(100);
        // print_bonds_including_H(100, 87600);
        // print_IMP_bonds(100);
        // print_atom_details_to_file();
        // print_dihedrals_without_H(100, 0);

        parmtop.close();
        return 0;
    }

    void AmberTopology::assign_hyperparameters()
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
        nconstraint_bonds_ = pointers_[12];
        nconstraint_angles_ = pointers_[13];
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

    void AmberTopology::create_molecules()
    {
        int start_index = 0;
        for (const int n_atoms_in_this_molecule : atoms_per_molecule_)
        {
            std::vector<int> atom_indices_for_this_molecule;
            for (int i = start_index; i < start_index+n_atoms_in_this_molecule; i++)
            {
                atom_indices_for_this_molecule.push_back(i);
            }
            Molecule mol = Molecule(atom_indices_for_this_molecule);
            Molecule_list_.push_back(mol);
            start_index += n_atoms_in_this_molecule;
        }
        for (Molecule& mol: Molecule_list_)
        {
            for (std::vector<int> atom_list = mol.get_atom_indices(); const int i : atom_list)
            {
                Atom& atom = atom_list_[i];
                mol.add_atom_name(atom.get_name());
                std::string residue_name = atom.get_residue_name();
                const unsigned long int residue_number = atom.get_residue_number();
                std::string residue_name_number = residue_name + std::to_string(residue_number);
                if (const bool check_if_exists = mol.check_if_residue_exists(residue_name_number); !check_if_exists)
                {
                    mol.add_residue_name(residue_name_number);
                }
            }
        }
    }

    int AmberTopology::read_coords(std::ifstream& coordfile) {
        std::string line;
        if (!std::getline(coordfile, line)) return 1;
        if (!std::getline(coordfile, line)) return 1; // Return files less than 2 lines in length.

        const std::size_t natoms = atom_list_.size();
        std::vector<double> box_dims;
        coordinates_.assign(atom_list_.size()*3, 0.0);

        int coord_processed = 0;
        while (std::getline(coordfile, line)) {
            if (util::check_if_line_empty(line)) continue;
            for (const std::vector<std::string> entries = util::split_line_over_empty_spaces(line); const std::string &entry: entries) {
                if (coord_processed < 3 * natoms) {
                    coordinates_[coord_processed] = (std::stod(entry));
                    coord_processed++;
                }
                else
                {
                    box_dims.push_back(std::stod(entry));
                }
            }
        }

        if (const size_t needed = natoms * 3; coord_processed < needed)
        {
            std::cerr << "Number of coordinates + box dimensions dont match the number of coordinates found in file. \n"
            << " Number of atoms in topology: " << natoms << "Number of coordinates found: " << (coordinates_.size())/3;
            return 1; // Return files that
        }

        if (box_dims.size() < 6)
        {
            std::cout << "Warning: Box dimensions not found in coordinate file. Assuming non-periodic.\n";
            box_x_ = 0.0; box_y_ = 0.0; box_z_ = 0.0;
            box_alpha_ = 90.0; box_beta_ = 90.0; box_gamma_ = 90.0; // 90 degrees is standard default
            return 0;
        }


        // Last six values are box dimensions
        box_x_ = box_dims[0];
        box_y_ = box_dims[1];
        box_z_ = box_dims[2];
        box_alpha_ = box_dims[3];
        box_beta_ = box_dims[4];
        box_gamma_ = box_dims[5];


        std::cout << "box: " << box_x_ << " " << box_y_ << " " << box_z_ << "\n" << box_alpha_ << " " << box_beta_ << " " << box_gamma_ << "\n";
        std::cout << "First three coordinates: " << coordinates_[0] << " " << coordinates_[1] << " " << coordinates_[2] << std::endl;
        return 0;
    }
} // namespace md