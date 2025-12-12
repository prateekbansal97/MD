
#include <iostream>
#include <string>
#include "include/AmberTopology/topology.h"
#include "include/System/System.h"
int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::cerr << "Usage: ./md <prmtop> [coords]" << std::endl;
        return 1;
    }

    std::string parmtop_path = argv[1];
    std::string coords_path = "";

    if (argc >= 3) {
        coords_path = argv[2];
    }

    Topology top = Topology::read_topology_coordinates(parmtop_path, coords_path);

    System sys = System(top);

    sys.init();
//    sys.minimize();
    return 0;
}
