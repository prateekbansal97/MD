#ifndef TOPOLOGY
#define TOPOLOGY_H

#include <vector>
#include "atom.h"

class Topology 

{
public:
    Topology(const std::vector<Atom>& atoms)
        : atom_list(atoms) {}
    std::vector<Atom> atom_list;

    std::vector<Atom> get_atoms() const {
        return atom_list;
    }
    std::vector<int> get_pointers();
};


#endif // TOPOLOGY_H