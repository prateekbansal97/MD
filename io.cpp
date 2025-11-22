
#include <string>
#include <vector>

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

std::vector<std::string> split_line_over_empty_spaces(std::string& line)
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
    return tokens;
}
