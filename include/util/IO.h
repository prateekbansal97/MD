
// #include <iostream>


#ifndef IO_H
#define IO_H
#include <string>
#include <vector>
#include <fstream>



namespace md::util
{
    bool check_first_char(const std::string& line, char target);
    bool check_if_line_starts_with_string(const std::string& line, const std::string& target);
    bool check_if_line_empty(const std::string& line);
    bool check_if_whitespace_in_string(const std::string& line);
    std::ifstream open_file(const std::string& path, std::ios::openmode mode = std::ios::in);

    std::vector<std::string> split_line_fixed_length(const std::string& line, size_t field_length);
    std::string remove_whitespaces_from_string(const std::string& str);
    std::string extract_header(const std::string& line);
    std::vector<std::string> split_line_on_delimiter(const std::string& line, const std::string& delimiter);
    std::vector<std::string> split_line_over_empty_spaces(const std::string& line, bool check_if_negative = false);

    std::vector<std::vector<std::string>> split_vector_into_chunks(const std::vector<std::string>& entries, const std::vector<std::size_t>& chunks);
}


#endif // IO_H