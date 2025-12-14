
#ifndef ATOM_H
#define ATOM_H

#include <string>
#include <utility>
#include <vector>

class Atom 
{
public:
    Atom(): partial_charge(0.0L) {}

    explicit Atom (std::string  atomname)
        : partial_charge(0.0f), name(std::move(atomname)) {}

    Atom(std::string& atomtype, bool isType)
    : partial_charge(0.0f),
    type(std::move(atomtype)) {}

    Atom(std::string& atomtype, const double charge): partial_charge(charge),
                                               type(std::move(atomtype)),
                                               residue_number(0),
                                               atomic_number(0),
                                               mass(0),
                                               atom_type_index(0),
                                               nExcluded_Atoms(0) {}

    void set_partial_charge(const double charge) { partial_charge = charge; }
    void set_type(const std::string& atomtype) { type = atomtype; }
    void set_residue_name(const std::string& resname) { residue_name = resname; }
    void set_residue_number(const int resnumber) { residue_number = resnumber; }
    void set_atomic_number(const int atomic_numbe) { this->atomic_number = atomic_numbe; }
    void set_mass(const double mas) { this->mass = mas; }
    void set_element(const std::string& elemen) { this->element = elemen; }
    void set_atom_type_index(const unsigned long int atom_type_inde) { this->atom_type_index = atom_type_inde; }
    void set_nExcluded_Atoms(const int n) { nExcluded_Atoms = n; }
    [[nodiscard]] unsigned long int get_atom_type_index() const { return atom_type_index; }

    [[nodiscard]] long double get_partial_charge() const { return partial_charge; }
    std::string& get_type() { return type; }
    std::string& get_name() { return name; }
    std::string& get_residue_name() { return residue_name; }
    [[nodiscard]] unsigned long int get_residue_number() const { return residue_number; }
    [[nodiscard]] int get_atomic_number() const { return atomic_number; }
    [[nodiscard]] int get_nExcluded_Atoms() const { return nExcluded_Atoms; }
    std::vector<int>& get_excluded_atoms() { return excluded_atoms; }
    double& get_mass() { return mass; }
    std::string& get_element() { return element; }
    [[nodiscard]] unsigned long int get_type_index() const { return atom_type_index; }

private:

    long double partial_charge;
    std::string type;
    std::string name;
    std::string residue_name;
    unsigned long int residue_number{};
    int atomic_number{};
    double mass{};
    std::string element;
    unsigned long int atom_type_index{};
    int nExcluded_Atoms{};
    std::vector<int> excluded_atoms;
    // unsigned long int nonbonded_parm_index;
};

#endif // ATOM_H
