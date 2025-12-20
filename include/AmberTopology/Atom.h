
#ifndef ATOM_H
#define ATOM_H

#include <string>
#include <utility>
#include <vector>

class Atom 
{
public:
    Atom(): partial_charge_(0.0L) {}

    explicit Atom (std::string  atomname)
        : partial_charge_(0.0f), name_(std::move(atomname)) {}

    Atom(std::string& atomtype_, bool isType)
    : partial_charge_(0.0f),
    atomtype_(std::move(atomtype_)) {}

    Atom(std::string& atomtype_, const double charge): partial_charge_(charge),
                                               atomtype_(std::move(atomtype_)) {}

    void set_partial_charge(const double charge) { partial_charge_ = charge; }
    void set_type(const std::string& type) { atomtype_ = type; }
    void set_residue_name(const std::string& resname) { residue_name_ = resname; }
    void set_residue_number(const int resnumber) { residue_number_ = resnumber; }
    void set_atomic_number(const int atomic_numbe) { this->atomic_number_ = atomic_numbe; }
    void set_mass(const double mas) { this->mass_ = mas; }
    void set_element(const std::string& elemen) { this->element_ = elemen; }
    void set_atom_type_index(const unsigned long int atom_type_inde) { this->atom_type_index_ = atom_type_inde; }
    void set_nExcluded_Atoms(const int n) { nExcluded_Atoms_ = n; }

    [[nodiscard]] const std::string& get_type() const { return atomtype_; };
    [[nodiscard]] const std::string& get_name() const { return name_; };
    [[nodiscard]] const std::string& get_residue_name() const { return residue_name_; }
    [[nodiscard]] unsigned long int get_residue_number() const { return residue_number_; }
    [[nodiscard]] int get_atomic_number() const { return atomic_number_; }
    [[nodiscard]] int get_nExcluded_Atoms() const { return nExcluded_Atoms_; }
    [[nodiscard]] unsigned long int get_atom_type_index() const { return atom_type_index_; }
    [[nodiscard]] long double get_partial_charge() const { return partial_charge_; }
    [[nodiscard]] std::vector<int>& get_excluded_atoms() { return excluded_atoms_; }
    [[nodiscard]] double get_mass() const { return mass_; }
    [[nodiscard]] const std::string& get_element() const { return element_; }

private:

    long double partial_charge_;
    std::string atomtype_;
    std::string name_;
    std::string residue_name_;
    unsigned long int residue_number_{};
    int atomic_number_{};
    double mass_{};
    std::string element_;
    unsigned long int atom_type_index_{};
    int nExcluded_Atoms_{};
    std::vector<int> excluded_atoms_;
    // unsigned long int nonbonded_parm_index;
};

#endif // ATOM_H
