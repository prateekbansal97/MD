
// #include <iostream>
// #include <fstream>

#ifndef IO_H
#define IO_H
#include <string>
#include <vector>



bool check_first_char(const std::string& line, char target);
bool check_if_line_starts_with_string(const std::string& line, const std::string& target);
bool check_if_line_empty(const std::string& line);

std::string extract_header(const std::string& line);
std::vector<std::string> split_line_on_delimiter(const std::string& line, std::string& delimiter);
std::vector<std::string> split_line_over_empty_spaces(std::string& line);
#endif // IO_H