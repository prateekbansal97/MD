#ifndef TOPOLOGY_H
#define TOPOLOGY_H

#include <vector>
#include "atom.h"
#include "io.h"

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
    
    void assign_hyperparameters();

private:
    std::vector<unsigned long int> pointers_;
    std::vector<Atom> atom_list;
    std::vector<std::string> residue_labels_;
    std::vector<float> bond_force_constants_;
    std::vector<float> bond_force_equils_;
    std::vector<float> angle_force_constants_;
    std::vector<float> angle_force_equils_;

    
    unsigned long int processed_atoms_index = 0;
    int index_processed = 0;
    int res_num = 0;
    int last_residue_num = 0;

    // nbmatrix of shape (nTypes_, nTypes_)
    std::vector<std::vector<unsigned long int>> nbmatrix_;

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