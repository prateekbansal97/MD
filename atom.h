
#include <string>

class Atom 
{
public:
    Atom(std::string atomtype)
        : type(atomtype), partial_charge(0.0f) {}

    Atom(std::string atomtype, float charge)
        : type(atomtype), partial_charge(charge) {}

        float partial_charge;
        std::string type;
};
