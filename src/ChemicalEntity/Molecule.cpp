#include "ChemicalEntity/Molecule.h"
#include <string>
#include <algorithm>
#include <iostream>

bool Molecule::check_if_residue_exists(const std::string& name)
{
    if (const auto it = std::ranges::find(residue_names_, name); it == residue_names_.end())
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
    for (std::string& name: residue_names_)
    {
        std::cout << name << "\n";
    }
}