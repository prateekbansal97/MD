
#include <string>
#include <cstdlib>
#include <vector>
#include <cctype>
#include <algorithm>
#include <numeric> 
#include <stdexcept>
#include <iostream>
#include <fstream>
#include "io.h"

bool check_first_char(const std::string& line, char target)
 {

    // Returns true if the first character of line matches target
    if (line.empty()) {
        return false;
    }
    return line[0] == target;
}

bool check_if_line_starts_with_string(const std::string& line, const std::string& target)
{
    // Returns true if the line starts with the target string
    if (line.size() < target.size()) {
        return false;
    }
    return line.compare(0, target.size(), target) == 0;
}


bool check_if_line_empty(const std::string& line)
{
    // Returns true if the line is empty or contains only whitespace
    return std::all_of(line.begin(), line.end(),
                       [](unsigned char ch) {
                           return std::isspace(static_cast<unsigned char>(ch)) != 0;
                       });
}


std::ifstream open_file(const std::string& path, std::ios::openmode mode)
{
    std::ifstream f(path, mode);
    if (!f) {
        std::cerr << "Error opening file: " << path << "\n";
        std::exit(1); // or exit(EXIT_FAILURE)
    }
    return f; // moved out
}

std::string extract_header(const std::string& line)
{
    // The header is going to be the string after the first '%' and after the first space
    // For example, the line "%FLAG POINTERS" would return "POINTERS"
    // the % will always be the first character of the line
    size_t first_space = line.find(' ');
    if (first_space == std::string::npos) {
        return "";
    }
    size_t end = line.find('\n');
    return line.substr(first_space + 1, end); // Exclude the '%' character
}

bool check_if_whitespace_in_string(const std::string& line)
{
    return std::any_of(line.begin(), line.end(),
                       [](unsigned char ch) {
                           return std::isspace(static_cast<unsigned char>(ch)) != 0;
                       });
}

std::vector<std::string> split_line_fixed_length(const std::string& line, size_t field_length)
{
    std::vector<std::string> tokens;
    for (size_t i = 0; i < line.length(); i += field_length) {
        tokens.push_back(line.substr(i, field_length));
    }
    return tokens;
}

std::string remove_whitespaces_from_string(const std::string& str)
{
    std::string result;
    std::copy_if(str.begin(), str.end(), std::back_inserter(result),
                 [](unsigned char ch) { return !std::isspace(ch); });
    return result;
}

std::vector<std::string> split_line_on_delimiter(const std::string& line, std::string& delimiter)
{
    std::vector<std::string> tokens;
    size_t start = 0;
    size_t end = line.find(delimiter);

    while (end != std::string::npos) {
        tokens.push_back(line.substr(start, end - start));
        start = end + delimiter.length();
        end = line.find(delimiter, start);
    }
    tokens.push_back(line.substr(start));

    return tokens;
}

std::vector<std::string> split_line_over_empty_spaces(const std::string& line, bool check_if_negative)
{
    // Split line over empty spaces of any length. 
    std::vector<std::string> tokens;
    size_t start = 0;
    size_t end = line.find_first_of(" \t", start);
    while (end != std::string::npos) {
        if (end != start) { // Avoid adding empty tokens
            tokens.push_back(line.substr(start, end - start));
        }
        start = line.find_first_not_of(" \t", end);
        end = line.find_first_of(" \t", start);
    }
    if (start != std::string::npos) {
        tokens.push_back(line.substr(start));
    }
    if (check_if_negative)
    {
        for (auto& element: tokens)
        {
            if (std::stod(element) < 0)
            {
                std::cout << "Negative element: " << element << std::endl;
            }
        }   
    }
    return tokens;
}

std::vector<std::vector<std::string>> split_vector_into_chunks(const std::vector<std::string>& entries, const std::vector<std::size_t>& chunks)
{
    // entries will be an array like ["0", "3", "4", "0", "3", "4", "0", "3", "4", "0"]
    // chunks will be an array like {3, 3, 4} for the entries array
    // First ensure that the sum of elements in chunks is the same as the length of entries
    // Then do the chunking and return the output


    const std::size_t total = std::accumulate(chunks.begin(), chunks.end(), std::size_t{0});
    if (total != entries.size()) {
        throw std::invalid_argument("sum(chunks) must equal entries.size()");
    }

    std::vector<std::vector<std::string>> result;
    result.reserve(chunks.size());

    std::size_t index = 0;
    for (std::size_t c : chunks) {
        result.emplace_back(entries.begin() + index, entries.begin() + index + c);
        index += c;
    }
    return result;
}


