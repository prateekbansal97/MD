#include "AmberTopology/molecule.h"
#include <string>
#include <algorithm>
#include <iostream>

bool Molecule::check_if_residue_exists(const std::string& name)
{
    if (const auto it = std::ranges::find(residue_names, name); it == residue_names.end())
    {
        return false;
    }
    else
    {
        return true;
    }
}

void Molecule::print_residue_names()
{
    for (std::string& name: residue_names)
    {
        std::cout << name << "\n";
    }
}