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
    std::vector<Atom> atom_list;

    std::vector<Atom> get_atoms() const {
        return atom_list;
    }

    int read_topology();
    std::vector<int> get_pointers();
    
private:
    std::vector<int> pointers_;
};


#endif // TOPOLOGY_H