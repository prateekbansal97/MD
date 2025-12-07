
#include <iostream>
#include <fstream>
#include <string>
#include "atom.h"
#include "io.h"
#include "topology.h"

int main (int argc, char* argv[])
{
    if (argc < 2) return 1;
    Topology topology;
    int p = topology.read_topology(argv[1]);
    return p;
}
