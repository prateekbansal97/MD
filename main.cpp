
#include <iostream>
#include <fstream>
#include <string>
#include "atom.h"
#include "io.h"
#include "topology.h"

int main ()
{ 
    Topology topology;
    int p = topology.read_topology();
    return p;
}