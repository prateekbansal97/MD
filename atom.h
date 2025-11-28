
#ifndef ATOM_H
#define ATOM_H

#include <string>

class Atom 
{
public:

    Atom (std::string atomname)
        : name(atomname), type(""), partial_charge(0.0f) {}
    Atom(std::string atomtype, bool isType)
        : type(atomtype), partial_charge(0.0f), name("") {}

    Atom(std::string atomtype, float charge)
        : type(atomtype), partial_charge(charge) {}


    long double partial_charge;
    std::string type;
    std::string name;
    std::string residue_name;
    unsigned long int residue_number;
    int atomic_number;
    float mass;
    std::string element;
    unsigned long int atom_type_index;
    int nExcluded_Atoms;
    // unsigned long int nonbonded_parm_index;
};

#endif // ATOM_H
