

#include <vector>
#include "atom.h"

class Topology 
{
public:
    Topology(const std::vector<Atom>& atoms)
        : atom_list(atoms) {}

    std::vector<Atom> atom_list;
};

std::vector<int> Topology::get_pointers() 
{
    std::vector<int> pointers;
    for (size_t i = 0; i < 10; ++i) {
        pointers.push_back(static_cast<int>(i));
    }
    return pointers;
}