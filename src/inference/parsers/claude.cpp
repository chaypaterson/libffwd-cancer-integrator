#include <iostream>
#include <vector>
#include <string>
#include <sstream>

std::vector<size_t> parseStringToVector(const std::string& str) {
    std::vector<size_t> result;
    std::stringstream ss(str);
    char c;
    size_t num;

    // Skip the opening bracket
    ss >> c;

    while (ss.good()) {
        std::string numStr;
        ss >> c;
        while (c != ',' && c != ']') {
            numStr += c;
            ss >> c;
        }

        result.push_back(std::stoul(numStr));

        if (c == ']') {
            break;
        }
    }

    return result;
}

int main() {
    std::string str = "[7, 48, 0, 13, 2]";
    std::vector<size_t> vect = parseStringToVector(str);

    for (size_t num : vect) {
        std::cout << num << " ";
    }
    std::cout << std::endl;

    return 0;
}

