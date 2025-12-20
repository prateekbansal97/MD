#ifndef TOPOLOGY_H
#define TOPOLOGY_H

#include "atom.h"
#include "harmonicUB.h"
#include "HarmonicImproper.h"
#include "HarmonicBond.h"
#include "HarmonicAngle.h"
#include "CosineDihedral.h"
#include "CMapGroup.h"
#include "molecule.h"

#include <vector>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <array>
#include <iostream>



class Topology 

{
public:
    Topology() = default;
    explicit Topology(const std::vector<Atom>& atoms)
        : atom_list_(atoms) {}

    // Implement a destructor
    

    [[nodiscard]] std::vector<Atom>& get_atoms() {
        return atom_list_;
    }

    [[nodiscard]] unsigned long int get_nAtoms() const {return nAtoms_;}

    static Topology read_topology_coordinates(const std::string& parmtop_path, const std::string& coords_path="");

    [[nodiscard]] size_t get_num_atoms() const { return atom_list_.size(); }

    std::vector<double>& get_coordinates() { return coordinates;}
    std::vector<HarmonicBond>& get_harmonic_bonds() {return HarmonicBond_list_;}
    std::vector<HarmonicAngle>& get_harmonic_angles() {return HarmonicAngle_list_;}
    std::vector<CosineDihedral>& get_cosine_dihedrals() {return CosineDihedral_list_;}
    std::vector<HarmonicUB>& get_harmonic_UBs() { return HarmonicUB_list_; }
    std::vector<HarmonicImproper>& get_harmonic_impropers() { return HarmonicIMP_list_; }
    std::vector<CMapGroup>& get_cmaps() { return CMapGroup_list_; }
    [[nodiscard]] const std::vector<int>& get_charmm_cmap_resolutions() const { return charmm_cmap_resolution_; }
    [[nodiscard]] const std::vector<double>& get_cmap_grid_data(const int set_index) const
    {
        switch (set_index) {
            case 1:
                return charmm_cmap_parameter_01;
            case 2:
                return charmm_cmap_parameter_02;
            case 3:
                return charmm_cmap_parameter_03;
            case 4:
                return charmm_cmap_parameter_04;
            case 5:
                return charmm_cmap_parameter_05;
            default:
                throw std::out_of_range("CMAP set index must be between 1 and 5. Received: " + std::to_string(set_index));
        }
    }

    std::vector<std::vector<unsigned long int>>& get_nb_matrix() { return nbmatrix_;}
    std::vector<double>& get_lennard_jones_Acoefs_() { return lennard_jones_Acoefs_; }
    std::vector<double>& get_lennard_jones_Bcoefs_() { return lennard_jones_Bcoefs_; }
    std::vector<double>& get_lennard_jones_14_Acoefs_() { return lennard_jones_14_Acoefs_; }
    std::vector<double>& get_lennard_jones_14_Bcoefs_() { return lennard_jones_14_Bcoefs_; }
    [[nodiscard]] std::array<double, 3> get_box_dimensions() const noexcept { return {box_x, box_y, box_z}; }

    int read_topology(std::ifstream& parmtop);

    [[maybe_unused]] std::vector<unsigned long int> get_pointers();

    [[maybe_unused]] [[maybe_unused]] void print_atom_details(int max_print);

    [[maybe_unused]] [[maybe_unused]] void print_atom_details_to_file();
    std::vector<Atom> atom_list_;
    void process_pointers_section(const std::string& line);
    void process_atom_names_section(const std::string& line);
    void process_charge_section(const std::string& line);
    void process_atomic_number_section(const std::string& line);
    void process_mass_section(const std::string& line);
    void process_atom_type_index_section(const std::string& line);
    void process_number_excluded_atoms_section(const std::string& line);
    void process_nonbonded_parm_index_section(const std::string& line);
    
    void process_residue_label_section(const std::string& line);
    void process_residue_pointer_section(const std::string& line);
    
    void process_bond_force_constant_section(const std::string& line);
    void process_bond_equil_section(const std::string& line);
    
    void process_angle_force_constant_section(const std::string& line);
    void process_angle_equil_section(const std::string& line);
   
    void process_charmm_urey_bradley_count_section(const std::string& line);
    void process_charmm_urey_bradley(const std::string& line);
    void process_charmm_urey_bradley_assign();
    void process_charmm_urey_bradley_force_constant(const std::string& line);
    void process_charmm_urey_bradley_equil_value(const std::string& line);
    void process_charmm_urey_bradley_assign_ffparams();
    [[maybe_unused]] void print_UB_bonds(int max_print);

    void process_dihedral_force_constant(const std::string& line);
    void process_dihedral_periodicity(const std::string& line);
    void process_dihedral_phase(const std::string& line);

    void process_scee_scale_factor(const std::string& line);
    void process_scnb_scale_factor(const std::string& line);

    void process_charmm_num_impropers(const std::string& line);
    void process_charmm_impropers(const std::string& line);
    void process_charmm_improper_assign();
    void process_charmm_num_impr_types(const std::string& line);
    void process_charmm_improper_force_constant(const std::string& line);
    void process_charmm_improper_phase(const std::string& line);
    void process_improper_bradley_assign_ffparams(); 
    [[maybe_unused]] void print_IMP_bonds(int max_print);

    void process_Lennard_Jones_Acoef(const std::string& line);
    void process_Lennard_Jones_Bcoef(const std::string& line);
    void process_Lennard_Jones_14_Acoef(const std::string& line);
    void process_Lennard_Jones_14_Bcoef(const std::string& line);

    void process_bonds_including_H(const std::string& line);// BONDS_INC_HYDROGEN
    void create_bonds_including_H();
    [[maybe_unused]] void print_bonds_without_H(int max_print, int start_point);
    void process_bonds_without_H(const std::string& line);
    void create_bonds_without_H();

    [[maybe_unused]] [[maybe_unused]] void print_bonds_including_H(int max_print, int start_point);

    void process_angles_including_H(const std::string& line);// ANGLES_INC_HYDROGEN
    void create_angles_including_H();

    [[maybe_unused]] [[maybe_unused]] void print_angles_without_H(int max_print, int start_point);
    void process_angles_without_H(const std::string& line);
    void create_angles_without_H();

    [[maybe_unused]] [[maybe_unused]] void print_angles_including_H(int max_print, int start_point);


    void process_dihedrals_including_H(const std::string& line);// DIHEDRALS_INC_HYDROGEN
    void create_dihedrals_including_H();
    [[maybe_unused]] void print_dihedrals_without_H(int max_print, int start_point);
    void process_dihedrals_without_H(const std::string& line);
    void create_dihedrals_without_H();
    [[maybe_unused]] void print_dihedrals_including_H(int max_print, int start_point);

    void process_excluded_atoms_list(const std::string& line);
    void create_excluded_atoms_list();

    void process_amber_atom_type(const std::string& line);
    
    void process_Charmm_Cmap_Count(const std::string& line);
    void process_Charmm_Cmap_Resolution(const std::string& line);
    void process_Charmm_Cmap_parameter_01(const std::string& line);
    void process_Charmm_Cmap_parameter_02(const std::string& line);
    void process_Charmm_Cmap_parameter_03(const std::string& line);
    void process_Charmm_Cmap_parameter_04(const std::string& line);
    void process_Charmm_Cmap_parameter_05(const std::string& line);
    void process_Charmm_Cmap_Index(const std::string& line);
    void create_Charmm_Cmap_Index_assign();
    void create_Cmap_Coefficient_Matrix_bicubic_spline();
    std::vector<double>& get_Cmap_Coefficient_Matrix_bicubic_spline() {return coeff_matrix_24x24;}

    void process_solvent_pointers(const std::string& line);
    
    void process_atoms_per_molecule(const std::string& line);
    void create_molecules();

    void process_box_dimensions(const std::string& line);

    void process_radius_set(const std::string& line);
    void process_radii(const std::string& line);
    void process_screen(const std::string& line);

    void process_ipol(const std::string& line);

    void build_lj14_pairlist();
    [[nodiscard]] bool is_14_pair(int i, int j) const;
    [[nodiscard]] double get_scee_scale_for_pair(int i, int j) const;


    int read_coords(std::ifstream& coordfile);
    template<typename T, typename... Args>
    void check_if_valid_indices(const std::string& vector_name, const std::vector<T>& vector, Args... args);

    void assign_hyperparameters();

    

private:
    std::vector<unsigned long int> pointers_;

    std::vector<std::string> residue_labels_;
    
    std::vector<double> bond_force_constants_;
    std::vector<double> bond_force_equils_;
    
    std::vector<double> angle_force_constants_;
    std::vector<double> angle_force_equils_;

    std::vector<int> charmm_urey_bradley_indices_;
    std::vector<HarmonicUB> HarmonicUB_list_;
    std::vector<double> charmm_urey_bradley_force_constants_;
    std::vector<double> charmm_urey_bradley_equil_values_;

    std::vector<double> dihedral_force_constants_;
    std::vector<double> dihedral_periodicities_;
    std::vector<double> dihedral_phases_;

    std::vector<double> scee_scale_factors_;
    std::vector<double> scnb_scale_factors_;

    std::vector<int> charmm_improper_indices_;
    std::vector<HarmonicImproper> HarmonicIMP_list_;
    std::vector<double> charmm_improper_force_constants_;
    std::vector<double> charmm_improper_phases_;

    std::vector<double> lennard_jones_Acoefs_;
    std::vector<double> lennard_jones_Bcoefs_;
    std::vector<double> lennard_jones_14_Acoefs_;
    std::vector<double> lennard_jones_14_Bcoefs_;

    std::vector<int> bonds_including_h_;
    std::vector<HarmonicBond> HarmonicBond_list_;
    std::vector<int> bonds_without_h_;

    std::vector<int> angles_including_h_;
    std::vector<HarmonicAngle> HarmonicAngle_list_;
    std::vector<int> angles_without_h_;

    std::vector<int> dihedrals_including_h_;
    std::vector<CosineDihedral> CosineDihedral_list_;
    std::vector<int> dihedrals_without_h_;

    std::vector<int> excluded_atoms_list_;

    std::vector<std::string> amber_atom_type_;
    
    std::vector<int> charmm_cmap_resolution_;
    std::vector<double> charmm_cmap_parameter_01; 
    std::vector<double> charmm_cmap_parameter_02;
    std::vector<double> charmm_cmap_parameter_03;
    std::vector<double> charmm_cmap_parameter_04;
    std::vector<double> charmm_cmap_parameter_05;
    std::vector<int> charmm_cmap_index;
    std::vector<double> coeff_matrix_24x24;
    std::vector<CMapGroup> CMapGroup_list_;

    std::vector<int> atoms_per_molecule_;
    std::vector<Molecule> Molecule_list;

    std::vector<double> radii;
    std::vector<double> screen;

    std::vector<double> coordinates;
    std::unordered_set<uint64_t> lj14_pairs_;
    std::unordered_map<uint64_t, double> scee14_scale_map_;


    unsigned long int processed_atoms_index = 0;
    int index_processed = 0;
    int res_num = 0;
    unsigned long int last_residue_num = 0;


    // nbmatrix of shape (nTypes_, nTypes_)
    std::vector<std::vector<unsigned long int>> nbmatrix_;

    unsigned long int nUreyBradley_{};
    unsigned long int nTypesUreyBradley{};

    unsigned long int nCharmmImpropers_{};
    unsigned long int nCharmmImproperTypes_{};

    // Hyperparameters
    unsigned long int nAtoms_{};
    unsigned long int nTypes_{};
    unsigned long int nBondHs_{};
    unsigned long int nBondAs_{};
    unsigned long int nAnglesH_{};
    unsigned long int nAnglesA_{};
    unsigned long int nTorsionsH_{};
    unsigned long int nTorsionsA_{};
    unsigned long int nexclusions_{};
    unsigned long int nresidues_{};
    unsigned long int nconstraint_bonds{};
    unsigned long int nconstraint_angles{};
    unsigned long int nconstraint_torsions_{};
    unsigned long int nunique_bond_types_{};
    unsigned long int nunique_angle_types_{};
    unsigned long int nunique_torsion_types_{};
    unsigned long int nsolty_terms_{};
    unsigned long int nphb_types_{};
    unsigned long int ifpert_{};
    unsigned long int nperturbed_bonds_{};
    unsigned long int nperturbed_angles_{};
    unsigned long int nperturbed_torsions_{};
    unsigned long int mperturbed_bonds_{};
    unsigned long int mperturbed_angles_{};
    unsigned long int mperturbed_torsions_{};
    unsigned long int ifbox_{};
    unsigned long int nmaxresidue_atoms_{};
    unsigned long int nCmap_{};
    unsigned long int nCmap_unique{};
    unsigned long int protein_res_{};
    unsigned long int nsolvent_mols_{};
    unsigned long int natoms_per_solvent_{};

    double box_beta{};
    double box_alpha{};
    double box_gamma{};
    double box_x{};
    double box_y{};
    double box_z{};

    std::string radius_set;

    int polarizable{};

};

template<typename T, typename... Args>
void Topology::check_if_valid_indices(const std::string& vector_name, const std::vector<T>& vector, Args... args)
{    
    auto check_index = [&](auto index) {
        if (static_cast<std::size_t>(index) >= vector.size() || static_cast<std::size_t>(index) < 0)
        {
            std::cerr << "Invalid index (" << index << ") passed. "
                     << "for vector " << vector_name 
                      << "Vector size is " << vector.size() << std::endl;
        }
    };
    
    ( (void)check_index(args), ... );
}



#endif // TOPOLOGY_H
