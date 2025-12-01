#ifndef TOPOLOGY_H
#define TOPOLOGY_H

#include <vector>
#include "atom.h"
#include "io.h"
#include "harmonicUB.h"
#include "HarmonicImproper.h"
#include "HarmonicBond.h"

class Topology 

{
public:
    Topology() {};
    Topology(const std::vector<Atom>& atoms)
        : atom_list(atoms) {}

    // Implement a destructor
    

    std::vector<Atom> get_atoms() const {
        return atom_list;
    }

    int read_topology();
    std::vector<unsigned long int> get_pointers();
    void print_atom_details(int max_print);
    void print_atom_details_to_file();
    
    void process_pointers_section(std::string& line);
    void process_atom_names_section(std::string& line);
    void process_charge_section(std::string& line);
    void process_atomic_number_section(std::string& line);
    void process_mass_section(std::string& line);
    void process_atom_type_index_section(std::string& line);
    void process_number_excluded_atoms_section(std::string& line);
    void process_nonbonded_parm_index_section(std::string& line);
    
    void process_residue_label_section(std::string& line);
    void process_residue_pointer_section(std::string& line);
    
    void process_bond_force_constant_section(std::string& line);
    void process_bond_equil_section(std::string& line);
    
    void process_angle_force_constant_section(std::string& line);
    void process_angle_equil_section(std::string& line);
   
    void process_charmm_urey_bradley_count_section(std::string& line);
    void process_charmm_urey_bradley(std::string& line);
    void process_charmm_urey_bradley_assign();
    void process_charmm_urey_bradley_force_constant(std::string& line);
    void process_charmm_urey_bradley_equil_value(std::string& line);
    void process_charmm_urey_bradley_assign_ffparams();
    void print_UB_bonds(int max_print);

    void process_dihedral_force_constant(std::string& line);
    void process_dihedral_periodicity(std::string& line);
    void process_dihedral_phase(std::string& line);

    void process_scee_scale_factor(std::string& line);
    void process_scnb_scale_factor(std::string& line);

    void process_charmm_num_impropers(std::string& line);
    void process_charmm_impropers(std::string& line);
    void process_charmm_improper_assign();
    void process_charmm_num_impr_types(std::string& line);
    void process_charmm_improper_force_constant(std::string& line);
    void process_charmm_improper_phase(std::string& line);
    void process_improper_bradley_assign_ffparams(); 
    void print_IMP_bonds(int max_print);

    void process_Lennard_Jones_Acoef(std::string& line);
    void process_Lennard_Jones_Bcoef(std::string& line);
    void process_Lennard_Jones_14_Acoef(std::string& line);
    void process_Lennard_Jones_14_Bcoef(std::string& line);

    void process_bonds_including_H(std::string& line);// BONDS_INC_HYDROGEN
    void process_bonds_including_H();
    void print_bonds(int max_print, int start_point);

    void assign_hyperparameters();

private:
    std::vector<unsigned long int> pointers_;
    std::vector<Atom> atom_list;
    std::vector<std::string> residue_labels_;
    
    std::vector<double> bond_force_constants_;
    std::vector<double> bond_force_equils_;
    
    std::vector<double> angle_force_constants_;
    std::vector<double> angle_force_equils_;

    std::vector<int> charmm_urey_bradley_indices;
    std::vector<HarmonicUB> HarmonicUB_list;
    std::vector<double> charmm_urey_bradley_force_constants_;
    std::vector<double> charmm_urey_bradley_equil_values_;

    std::vector<double> dihedral_force_constants_;
    std::vector<double> dihedral_equil_values_;
    std::vector<double> dihedral_periodicities_;
    std::vector<double> dihedral_phases_;

    std::vector<double> scee_scale_factors_;
    std::vector<double> scnb_scale_factors_;

    std::vector<int> charmm_improper_indices_;
    std::vector<HarmonicImproper> HarmonicIMP_list;
    std::vector<double> charmm_improper_force_constants_;
    std::vector<double> charmm_improper_phases_;

    std::vector<double> lennard_jones_Acoefs_;
    std::vector<double> lennard_jones_Bcoefs_;
    std::vector<double> lennard_jones_14_Acoefs_;
    std::vector<double> lennard_jones_14_Bcoefs_;

    std::vector<int> bonds_including_h_;
    std::vector<HarmonicBond> HarmonicBond_list;

    unsigned long int processed_atoms_index = 0;
    int index_processed = 0;
    int res_num = 0;
    int last_residue_num = 0;

    // nbmatrix of shape (nTypes_, nTypes_)
    std::vector<std::vector<unsigned long int>> nbmatrix_;

    unsigned long int nUreyBradley_;
    unsigned long int nTypesUreyBradley;

    unsigned long int nCharmmImpropers_;
    unsigned long int nCharmmImproperTypes_;
    // Hyperparameters
    unsigned long int nAtoms_;
    unsigned long int nTypes_;
    unsigned long int nBondHs_;
    unsigned long int nBondAs_;
    unsigned long int nAnglesH_;
    unsigned long int nAnglesA_;
    unsigned long int nTorsionsH_;
    unsigned long int nTorsionsA_;
    unsigned long int nexclusions_;
    unsigned long int nresidues_;
    unsigned long int nconstraint_bonds;
    unsigned long int nconstraint_angles;
    unsigned long int nconstraint_torsions_;
    unsigned long int nunique_bond_types_;
    unsigned long int nunique_angle_types_;
    unsigned long int nunique_torsion_types_;
    unsigned long int nsolty_terms_;
    unsigned long int nphb_types_;
    unsigned long int ifpert_;
    unsigned long int nperturbed_bonds_;
    unsigned long int nperturbed_angles_;
    unsigned long int nperturbed_torsions_;
    unsigned long int mperturbed_bonds_;
    unsigned long int mperturbed_angles_;
    unsigned long int mperturbed_torsions_;
    unsigned long int ifbox_;
    unsigned long int nmaxresidue_atoms_;

};


#endif // TOPOLOGY_H