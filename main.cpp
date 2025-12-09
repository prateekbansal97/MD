
#include <iostream>
#include <fstream>
#include <string>
#include "include/AmberTopology/atom.h"
#include "include/AmberTopology/io.h"
#include "include/AmberTopology/topology.h"

int main(int argc, char* argv[])
{
    if (argc < 2) return 1;
    Topology topology;

    if (argc >= 3) return topology.read_topology_coordinates(argv[1], argv[2]);
    return topology.read_topology_coordinates(argv[1]);
}
