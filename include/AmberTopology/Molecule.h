
#ifndef MOLECULE_H
#define MOLECULE_H

#include <string>
#include <vector>

class Molecule 
{
public:
    explicit Molecule(const std::vector<int>& atom_indices):
    atom_indices_(atom_indices),
    nResidues_(0)
    {
        nAtoms_ = atom_indices.size();
    }

    void set_molecule_name(const std::string& name) {this->molecule_name_ = name;}
    void add_atom_name(const std::string& name) {atom_names_.push_back(name);}
    void add_residue_name(const std::string& name) {residue_names_.push_back(name);}

    [[nodiscard]] const std::string& get_molecule_name() {return molecule_name_;}
    [[nodiscard]] const std::vector<int>& get_atom_indices() {return atom_indices_;}
    [[nodiscard]] bool check_if_residue_exists(const std::string& name);

    void print_residue_names();

private:
    std::vector<int> atom_indices_;
    std::vector<std::string> atom_names_;
    std::vector<std::string> residue_names_;
    std::string molecule_name_;
    unsigned long int nAtoms_;
    unsigned long int nResidues_;
};

#endif