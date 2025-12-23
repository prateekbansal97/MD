
#include <iostream>
#include <string>
#include "../include/AmberTopology/AmberTopology.h"
#include "../include/System/System.h"

using namespace md::AmberTopology;

int main(const int argc, char* argv[])
{
    if (argc < 2) {
        std::cerr << "Usage: ./md <prmtop> [coords]" << std::endl;
        return 1;
    }

    const std::string parmtop_path = argv[1];
    std::string coords_path;

    if (argc >= 3) {
        coords_path = argv[2];
    }

    const AmberTopology top = AmberTopology::read_topology_coordinates(parmtop_path, coords_path);

    md::System sys = md::System(top);

    sys.init();
//    sys.minimize();
    return 0;
}
