
#ifndef MOLECULE_H
#define MOLECULE_H

#include <string>
#include <vector>

class Molecule 
{
public:
    Molecule(std::vector<int> atom_indices):
    atom_indices(atom_indices),
    molecule_name(""),
    nResidues(0) 
    {
        nAtoms = atom_indices.size();
    }

    void set_molecule_name(std::string& name) {this->molecule_name = name;}
    std::string get_molecule_name() {return molecule_name;}

    void add_atom_name(std::string& name) {atom_names.push_back(name);}
    void add_residue_name(std::string& name) {residue_names.push_back(name);}

    
    std::vector<int> get_atom_indices() {return atom_indices;}

    bool check_if_residue_exists(std::string& name);
    void print_residue_names();

private:
    std::vector<int> atom_indices;
    std::vector<std::string> atom_names;
    std::vector<std::string> residue_names;
    std::string molecule_name;
    unsigned long int nAtoms;
    unsigned long int nResidues;    
};

#endif