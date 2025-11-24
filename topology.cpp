#include <iostream>
#include <fstream>
#include <string>
#include "topology.h"
#include <vector>
#include "atom.h"
#include "io.h"

std::vector<double> Topology::get_pointers() 
{
    return pointers_;
}

int Topology::read_topology()
{
    std::string filename = "/Users/Prateek/Desktop/Coding/MD/CB1_apo_assym.prmtop";
    std::ifstream file(filename);

    if (!file) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return 1;
    }

    double line_number = 0;
    std::string line;
    while (std::getline(file, line)) {

        bool pointers_processed = false;
        bool atom_names_processed = false;
        // Look specifically for the POINTERS flag
        if (check_first_char(line, '%'))
        {
            std::string header = extract_header(line);
            std::cout << "Header: " << header << std::endl;
            if (header == "POINTERS")
            {
                bool pointers_processed = true;
            }
            else if (header == "ATOM_NAME")
            {
                bool atom_names_processed = true;
            }
            else {
                std::cout << "Skipping unhandled header: " << header << std::endl;
                continue;
            }
        }
        else if (check_if_line_starts_with_string(line, "%FORMAT"))
        {
            continue;
        }
        else if (check_if_line_empty(line))
        {
            continue;
        }
        else {
            if (pointers_processed)
            {
                process_pointers_section(line);
            }
            else if (atom_names_processed)
            {
                process_atom_name_section(line);
            }
        }
    }

    file.close();
    return 0;
}
void Topology::process_pointers_section(std::string& line)
{
    std::vector<std::string> entries = split_line_over_empty_spaces(line);
    for (const auto& entry : entries) {
        try {
            double value = std::stod(entry);
            pointers_.push_back(value);
        } catch (const std::invalid_argument& e) {
            std::cerr << "Invalid number format: " << entry << std::endl;
        } catch (const std::out_of_range& e) {
            std::cerr << "Number out of range: " << entry << std::endl;
        }
}

void Topology::process_atom_name_section(std::string& line)
{
    std::vector<std::string> entries = split_line_over_empty_spaces(line);
    for (const auto& entry : entries) {
        Atom atom(entry);
        atom_list.push_back(atom);
    }
}
