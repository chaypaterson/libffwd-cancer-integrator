#include <iostream>
#include <vector>
#include <sstream>
#include <string>

std::vector<size_t> parseStringToVector(const std::string& str) {
    std::vector<size_t> result;
    std::stringstream ss(str);
    char c;
    size_t num;

    ss >> c; // consume the opening '['
    while (ss >> num) {
        result.push_back(num);
        ss >> c; // consume the comma or closing ']'
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
