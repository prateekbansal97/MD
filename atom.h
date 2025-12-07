
#ifndef ATOM_H
#define ATOM_H

#include <string>
#include <vector>

class Atom 
{
public:
    Atom(): partial_charge(0.0L), type(""), name(""), residue_name(""), 
            residue_number(0), atomic_number(0), 
            mass(0.0f), element(""), 
            atom_type_index(0), nExcluded_Atoms(0) {}
    
    Atom (std::string atomname)
        : name(atomname), type(""), partial_charge(0.0f) {}
    Atom(std::string atomtype, bool isType)
        : type(atomtype), partial_charge(0.0f), name("") {}

    Atom(std::string atomtype, double charge)
        : type(atomtype), partial_charge(charge) {}


    long double partial_charge;
    std::string type;
    std::string name;
    std::string residue_name;
    unsigned long int residue_number;
    int atomic_number;
    double mass;
    std::string element;
    unsigned long int atom_type_index;
    int nExcluded_Atoms;
    std::vector<int> excluded_atoms;
    // unsigned long int nonbonded_parm_index;
};

#endif // ATOM_H
