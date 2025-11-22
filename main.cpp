
#include <iostream>
#include <fstream>
#include <string>
#include "atom.h"
#include "io.h"
// #include "topology.h"

int main ()
{ 
    std::string filename = "/Users/Prateek/Desktop/Coding/MD/CB1_apo_assym.prmtop";
    // std::cout << "Enter the filename: ";
    // std::cin >> filename;

    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return 1;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::string target = "%FLAG POINTERS";
        if (check_if_line_starts_with_string(line, target)) {
            std::cout << line << std::endl;

            std::string dataLine;
            std::string tmpLine;
            int steps = 2; // skip the next line and then advance 12 more lines
            bool ok = true;
            for (int i = 0; i < steps; ++i) {
                if (!std::getline(file, tmpLine)) { ok = false; break; }
            }
            if (ok) {
                dataLine = tmpLine;
            } else {
                std::cerr << "Reached EOF before finding the requested line\n";
            }
            std::cout << "Data Line: " << dataLine << std::endl;
    }
}
    file.close();
    return 0;
}